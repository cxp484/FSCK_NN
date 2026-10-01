! Standalone FSK V4 binary-table exporter. No OpenFOAM, ONNX, or Python required.
! Usage: FSCK_to_csv DATABASE_ROOT P_MIN_ATM P_MAX_ATM OUTPUT.csv
! Exports existing files on the declared V4 composition/pressure grids.
!
! Build from FSCK_NN:
!   gfortran -std=f2008 -O2 test/FSCK_to_csv.f90 -o test/build/FSCK_to_csv
! Example (export all existing V4 files from 1 through 5 atm):
!   test/build/FSCK_to_csv /path/to/FSKTableV4 1 5 output.csv
!
! Input format, for each fixed pressure/composition file:
!   28 local temperatures x 28 reference temperatures = 784 records.
!   Each record: 32 quadrature kappas + kappa_max, all native-endian float32.
!   Tg is the outer index; Tr is the inner (fastest-varying) index.
! Output: one CSV row per temperature pair, with 39 columns:
!   P, xCO2, xH2O, xCO, fv, Tg, Tr, kappa_01, ..., kappa_32.
! No interpolation, NN evaluation, or unit conversion is performed.
program FSCK_to_csv
  use iso_fortran_env, only: real32, real64, int64, error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  ! These axes match PDB, xCO2DB, xH2ODB, xCODB, and fvDB in fskTableV4.f90.
  ! The exporter constructs known V4 filenames rather than listing directories.
  ! Missing combinations are skipped; custom axes require updating these arrays.
  ! Units: pressure in atm; gases as mole fractions; soot as volume fraction.
  real(real64), parameter :: pressures(*) = [ &
    .1d0,.2d0,.3d0,.4d0,.5d0,.7d0,1d0,2d0,3d0,4d0,5d0,6d0,7d0, &
    8d0,9d0,10d0,11d0,12d0,13d0,14d0,15d0,20d0,25d0,30d0,35d0, &
    40d0,45d0,50d0,55d0,60d0,65d0,70d0,75d0,80d0]
  real(real64), parameter :: co2(*) = [0d0,.01d0,.02d0,.03d0,.04d0,.05d0,.25d0,.5d0,.75d0,1d0]
  real(real64), parameter :: h2o(*) = [0d0,.01d0,.02d0,.03d0,.04d0,.05d0,.1d0,.15d0,.2d0,.25d0,.5d0,.75d0,1d0]
  real(real64), parameter :: co(*) = [0d0,.01d0,.05d0,.1d0,.25d0,.5d0]
  real(real64), parameter :: soot(*) = [0d0,1d-7,1d-6,1d-5]

  ! Deferred-length strings preserve complete command-line arguments and paths.
  ! suffix holds the fixed V4 filename; message holds a runtime I/O diagnostic.
  character(:), allocatable :: database,output,arg,filename
  character(256) :: suffix
  character(1024) :: message
  ! values is one CSV row: 1:5 = pressure/composition, 6:7 = temperatures,
  ! 8:39 = the 32 quadrature kappas. raw includes the unused final kappa_max.
  ! Memory use does not grow with the number of files or exported rows.
  real(real64) :: pmin,pmax,values(39)
  real(real32) :: raw(33)
  ! Loop indices: pressure, CO2, H2O, CO, soot, Tg, Tr, quadrature, character.
  ! NEWUNIT assigns free file-unit numbers; IOSTAT=0 indicates successful I/O.
  integer :: ip,ic,ih,ico,iv,it,ir,q,j,input_unit,output_unit,ios
  ! Use 64-bit counters because a full database export may contain many rows.
  integer(int64) :: bytes,nfiles,nrows,pressure_files
  ! Track opened files so fail() can close input and delete incomplete output.
  logical :: exists,output_open=.false.,input_open=.false.

  ! 1. Read and validate the requested inclusive pressure interval.
  ! Equal bounds select a single tabulated pressure. Intermediate bounds select
  ! the database pressures inside the interval, not interpolated pressure states.
  if (command_argument_count()/=4) then
    write(error_unit,'(a)') 'Usage: FSCK_to_csv DATABASE_ROOT P_MIN_ATM P_MAX_ATM OUTPUT.csv'
    write(error_unit,'(a)') 'Inclusive pressure range; use equal limits for one pressure. Output must not exist.'
    stop 1
  endif
  call argument(1,database)
  call argument(2,arg)
  read(arg,*,iostat=ios) pmin
  if (ios/=0) call fail('Invalid minimum pressure')
  call argument(3,arg)
  read(arg,*,iostat=ios) pmax
  if (ios/=0) call fail('Invalid maximum pressure')
  call argument(4,output)
  if (.not.ieee_is_finite(pmin).or..not.ieee_is_finite(pmax)) call fail('Pressures must be finite')
  if (pmin<=0d0.or.pmax<pmin) call fail('Require 0 < P_MIN <= P_MAX')
  if (.not.any(pressures>=pmin.and.pressures<=pmax)) call fail('No V4 pressures in the requested interval')
  if (storage_size(raw(1))/=32) call fail('A 32-bit real kind is required')

  ! 2. Create the CSV. STATUS='new' refuses to overwrite an existing file.
  ! A failed OPEN leaves output_open false, so fail() cannot delete that file.
  open(newunit=output_unit,file=output,status='new',action='write',iostat=ios,iomsg=message)
  if (ios/=0) call fail('Cannot create output: '//trim(message))
  output_open=.true.
  ! Nonadvancing writes assemble the header on one line; I2.2 yields 01...32.
  write(output_unit,'(a)',advance='no',iostat=ios) 'P,xCO2,xH2O,xCO,fv,Tg,Tr'
  if (ios/=0) call fail('Cannot write CSV header')
  do q=1,32
    write(output_unit,'(a,i2.2)',advance='no',iostat=ios) ',kappa_',q
    if (ios/=0) call fail('Cannot write CSV header')
  enddo
  write(output_unit,*,iostat=ios)
  if (ios/=0) call fail('Cannot write CSV header')

  ! 3. Visit every pressure/composition combination on the declared V4 axes.
  ! Export stored combinations as-is, without filtering by the sum of gas
  ! mole fractions. Pressure and composition are fixed within each binary file.
  nfiles=0; nrows=0
  do ip=1,size(pressures)
    if (pressures(ip)<pmin.or.pressures(ip)>pmax) cycle
    pressure_files=0
    do ic=1,size(co2)
    do ih=1,size(h2o)
    do ico=1,size(co)
    do iv=1,size(soot)
      ! Label 100 below reproduces loadIntoMem's filename format. Pressure
      ! occurs twice: once in the directory name and once in the filename.
      write(suffix,100) pressures(ip),pressures(ip),co2(ic),h2o(ih),co(ico),soot(iv)
      ! Match V4 loadIntoMem's strpad: replace field-padding spaces by zeros.
      ! Example suffix: /01.0/P.01.0_CO2.0.00_H2O.0.00_CO.0.00_fv.1.E-05_k.dat
      do j=1,len_trim(suffix)
        if (suffix(j:j)==' ') suffix(j:j)='0'
      enddo
      filename=database//trim(suffix)
      inquire(file=filename,exist=exists)
      if (.not.exists) cycle
      ! Expected size: 28 * 28 * 33 * 4 = 103488 bytes. The V4 format has
      ! neither a header nor sequential-unformatted record markers.
      inquire(file=filename,size=bytes,iostat=ios)
      if (ios/=0) call fail('Cannot inspect '//filename)
      if (bytes/=103488_int64) call fail('Expected 103488-byte V4 file: '//filename)
      ! Stream access reads the same contiguous bytes as the original direct-
      ! access records. No endian conversion is requested: database byte order
      ! must match this machine's native representation.
      open(newunit=input_unit,file=filename,access='stream',form='unformatted', &
           status='old',action='read',iostat=ios,iomsg=message)
      if (ios/=0) call fail('Cannot read '//filename//': '//trim(message))
      input_open=.true.
      values(1:5)=[pressures(ip),co2(ic),h2o(ih),co(ico),soot(iv)]
      ! 4. Read one 33-float record per temperature pair. The record number
      ! is 28*(it-1)+ir; both temperatures span 300 to 3000 K in 100 K steps.
      do it=1,28
      do ir=1,28
        read(input_unit,iostat=ios) raw
        if (ios/=0) call fail('Cannot read temperature record in '//filename)
        if (.not.all(ieee_is_finite(raw)).or.any(raw<0)) call fail('Invalid kappa data in '//filename)
        values(6)=300d0+100d0*(it-1)
        values(7)=300d0+100d0*(ir-1)
        ! Omit raw(33), which is kappa_max, not a quadrature value.
        ! Widen float32 to float64 before formatting, without adding precision
        ! to the source data. The CSV digits preserve the original float32 value.
        values(8:39)=real(raw(1:32),real64)
        ! Write one number followed by 38 comma-prefixed numbers. ES uses
        ! scientific notation with enough digits to preserve stored kappas.
        write(output_unit,'(es24.16e3,38(",",es24.16e3))',iostat=ios) values
        if (ios/=0) call fail('Cannot write CSV data (check available disk space)')
        nrows=nrows+1
      enddo
      enddo
      close(input_unit)
      input_open=.false.
      nfiles=nfiles+1; pressure_files=pressure_files+1
    enddo
    enddo
    enddo
    enddo
    write(*,'(a,f6.1,a,i0)') 'Pressure ',pressures(ip),' atm: files = ',pressure_files
  enddo

  ! 5. Reject an empty export, flush buffered data, and report completion.
  ! In particular, a wrong database path must not produce a successful-looking
  ! CSV containing only a header.
  if (nfiles==0) call fail('No matching database files found; check the database path and pressure limits')
  flush(output_unit,iostat=ios)
  if (ios/=0) call fail('Cannot flush CSV output')
  close(output_unit,iostat=ios)
  output_open=.false.
  if (ios/=0) call fail('Cannot close CSV output')
  write(*,'(a,i0,a,i0,a)') 'Exported ',nfiles,' files, ',nrows,' rows.'
  write(*,'(a)') 'CSV: '//output
  ! Fixed-width V4 filename fields; the zero-padding pass is applied above.
100 format('/',f4.1,'/P.',f4.1,'_CO2.',f4.2,'_H2O.',f4.2,'_CO.',f4.2,'_fv.',es6.0,'_k.dat')
contains
  ! Read an argument in two passes: first determine its length, then allocate
  ! exactly that many characters. Quote paths containing spaces in the shell.
  subroutine argument(index,text)
    integer, intent(in) :: index
    character(:), allocatable, intent(out) :: text
    integer :: length,status
    call get_command_argument(index,length=length,status=status)
    if (status/=0.or.length==0) call fail('Empty or invalid command argument')
    allocate(character(length)::text)
    call get_command_argument(index,text,status=status)
    if (status/=0) call fail('Cannot read command argument')
  end subroutine

  ! Centralized error exit. Report to standard error and return a nonzero
  ! process status. Delete output only while our newly created CSV is open;
  ! this prevents a partially written export from being mistaken for complete.
  subroutine fail(text)
    character(*), intent(in) :: text
    integer :: close_status
    write(error_unit,'(a)') 'FSCK_to_csv: '//text
    if (input_open) close(input_unit,iostat=close_status)
    if (output_open) close(output_unit,status='delete',iostat=close_status)
    stop 1
  end subroutine
end program FSCK_to_csv
