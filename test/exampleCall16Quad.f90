! Runtime example: request 16 quadrature values from the 32-point FSK V4 table.
!
! Build from FSCK_NN: bash test/build_exampleCall16Quad.sh
! Run: test/build/exampleCall16Quad /path/to/FSKTableV4
!
! This demonstrates the Fortran calls behind spectralModels.f90's
! fsk_gw(...,'V4Table') and fsk_v4_table_nq wrappers. It does not need OpenFOAM
! or ONNX. The example prints optical properties, not a complete RTE solution.
program exampleCall16Quad
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use commonData, only: GasMixInfo
  use FSKDBPath, only: v4path
  use quadlib, only: quadgen2
  use fskTableV4, only: loadIntoMem, fsk_v4_table, freeMem
  implicit none

  ! NoPDB is fixed by this database's record layout. Nq is the number of
  ! spectral transport equations/quadrature points chosen by the solver.
  integer, parameter :: NoPDB=32, Nq=16, nCells=2
  real(real64) :: gFSK(NoPDB),wFSK(NoPDB),g(Nq),w(Nq),k(Nq),a(Nq)
  real(real64) :: pressure_limits(2),gas_max(3),temperature_max,soot_max
  type(GasMixInfo) :: cells(nCells),reference
  character(1024) :: database
  integer :: cell,q,status,path_length

  if (command_argument_count()/=1) then
    print *, 'Usage: exampleCall16Quad DATABASE_ROOT'
    stop 1
  endif
  call get_command_argument(1,database,length=path_length,status=status)
  ! The original loader uses a 250-character buffer for the full filename.
  if (status/=0.or.path_length==0.or.path_length>150) error stop 'Invalid or excessively long database path'

  ! 1. Example cell states. A radiation solver obtains these from its fields.
  ! Pressure is in atm, temperatures in K, gases are mole fractions, and fv
  ! is soot volume fraction. These off-grid values exercise interpolation.
  cells(1)%P=1d0
  cells(1)%T=950d0
  cells(1)%xCO2=.015d0
  cells(1)%xH2O=.015d0
  cells(1)%xCO=.005d0
  cells(1)%fv=5d-8
  cells(1)%xCH4=0d0
  cells(1)%xC2H4=0d0
  cells(2)=cells(1)
  cells(2)%T=1450d0
  cells(2)%xH2O=.035d0

  ! A common reference temperature supplies the Planck-weighting convention.
  ! get_ka, called internally below, uses reference%T and the local mixture.
  reference=cells(1)
  reference%T=1100d0

  ! 2. Initialization: generate BOTH grids using the V4 quadrature convention.
  ! These calls are equivalent to fsk_gw(32,...,'V4Table') and
  ! fsk_gw(16,...,'V4Table'). g contains abscissae, w the integration weights.
  ! The 16 abscissae are not obtained by selecting alternate 32-point entries.
  call quadgen2(.false.,gFSK,wFSK,NoPDB,2d0)
  call quadgen2(.false.,g,w,Nq,2d0)

  ! Load the necessary database region ONCE, outside the cell/quadrature loop.
  ! In a transient solver these bounds should cover the states to be queried;
  ! reload if the solution leaves the loaded region. The original loader reads
  ! pressure brackets and composition prefixes, so bounding files must exist.
  v4path=trim(database)
  pressure_limits=[minval(cells%P),maxval(cells%P)]
  gas_max=[maxval(cells%xCO2),maxval(cells%xH2O),maxval(cells%xCO)]
  temperature_max=max(maxval(cells%T),reference%T)
  soot_max=maxval(cells%fv)
  call loadIntoMem(pressure_limits,temperature_max,gas_max,soot_max)

  ! 3. Runtime cell loop: the table is now in memory; this call performs
  ! interpolation and does not reread the binary files for each cell.
  do cell=1,nCells
    k=0d0; a=0d0
    call fsk_v4_table(cells(cell),reference,Nq,k,a,g,gFSK,wFSK,NoPDB)
    ! Internally fsk_v4_table:
    !   (a) get_ka computes kpre(32) and apre(32) at the cell state;
    !       a uses distributions weighted at reference%T and cells(cell)%T.
    !   (b) k_interp maps kpre from gFSK(32) to the requested g(16).
    !   (c) a_interp maps apre from gFSK(32) to the requested g(16).
    ! Thus a is formed on the database grid BEFORE mapping to 16 points.
    if (.not.all(ieee_is_finite(k)).or..not.all(ieee_is_finite(a))) error stop 'Nonfinite optical properties'

    ! Match the final kappa floor in spectralModels.f90:fsk_v4_table_nq.
    ! This floor is applied by the solver wrapper, not by fsk_v4_table itself.
    k=max(k,1d-9)
    write(*,'(/,a,i0,a,f8.1,a,f8.1)') 'Cell ',cell,': Tg = ',cells(cell)%T,', Tr = ',reference%T
    write(*,'(a)') ' q        g(q)             w(q)             kappa(q)           a(q)'
    do q=1,Nq
      write(*,'(i3,4(1x,es18.10))') q,g(q),w(q),k(q),a(q)
      ! A nonscattering thermal RTE solver uses k(q) for attenuation and
      ! k(q)*a(q)*Ib(Tg) as its emission source at this quadrature point.
      ! After solving the 16 RTEs, integrate intensity with sum(w(:)*I(:)).
      ! Use the 16-point weights w, not the 32-point database weights wFSK.
    enddo
    write(*,'(a,es18.10)') 'Sum of solver weights: ',sum(w)
    write(*,'(a,es18.10)') 'Discrete sum w*a:      ',sum(w*a)
    ! Finite quadrature, the original afun endpoint treatment, and the mapping
    ! can make sum(w*a) differ from one. No extra normalization is applied here.
  enddo

  ! 4. Finalization: release database memory when the solver is finished.
  call freeMem()
end program exampleCall16Quad
