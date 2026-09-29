! Added standalone entry point. Kappa-only path; original interpolation routines are unchanged.
subroutine fsck_file(path, pathlen, state, results, nodes, weights, status) bind(C, name='fsck_file')
  use iso_c_binding
  use commonData, only: GasMixInfo
  use quadlib, only: quadgen2
  use fskTable, only: k_interp
  implicit none
  integer(c_int), value :: pathlen
  character(c_char), intent(in) :: path(pathlen)
  real(c_double), intent(in) :: state(5) ! P[atm], xCO2, xH2O, xCO, fv
  real(c_double), intent(out) :: results(32,784), nodes(32), weights(32)
  integer(c_int), intent(out) :: status
  character(:), allocatable :: filename
  integer :: unitno, ios, j, local_i, ref_i, record
  real(4) :: raw(33)
  real(8) :: k(32), kpre(32)
  type(GasMixInfo) :: mix
  status=0
  results=0
  allocate(character(pathlen) :: filename)
  do j=1,pathlen
    filename(j:j)=path(j)
  enddo
  if (.not. allocated(kq)) then
    nP=1; nfv=1; nxco2=1; nxh2o=1; nxco=1; nTg=28; nT=28
    allocate(P(1),fv(1),xCO2(1),xH2O(1),xCO(1),Tg(28),T(28))
    allocate(kq(1,1,1,1,1,28,28))
    do j=1,28
      Tg(j)=300.d0+(j-1)*100.d0
      T(j)=Tg(j)
    enddo
  endif
  P(1)=state(1); xCO2(1)=state(2); xH2O(1)=state(3); xCO(1)=state(4); fv(1)=state(5)
  open(newunit=unitno,file=filename,access='direct',form='unformatted',recl=132, &
       status='old',action='read',iostat=ios)
  if (ios/=0) then
    status=1
    return
  endif
  do local_i=1,28
    do ref_i=1,28
      record=(local_i-1)*28+ref_i
      read(unitno,rec=record,iostat=ios) raw
      if (ios/=0) then
        close(unitno)
        status=2
        return
      endif
      ! Same selection as loadIntoMem for NDB=Nq=32: first 32 of 33 floats.
      kq(1,1,1,1,1,local_i,ref_i)%values=raw(1:32)
    enddo
  enddo
  close(unitno)
  ! Match spectralModels.f90 fsk_gw(...,"V4Table").
  call quadgen2(.false.,nodes,weights,32,2.d0)
  mix%P=state(1); mix%xCO2=state(2); mix%xH2O=state(3); mix%xCO=state(4)
  mix%fv=state(5); mix%xCH4=0; mix%xC2H4=0
  do local_i=1,28
    mix%T=Tg(local_i)
    do ref_i=1,28
      ! Same kappa path as fsk_v4_table, without computing the a-function.
      call get_k(mix,T(ref_i),kpre)
      call k_interp(32,nodes,kpre,32,nodes,k)
      record=(local_i-1)*28+ref_i
      results(1:32,record)=k
    enddo
  enddo
end subroutine fsck_file
