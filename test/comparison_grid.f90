module comparison_grid
  use iso_fortran_env, only: real64, int64, error_unit
  implicit none
  private
  public :: read_grid
contains
  subroutine read_grid(filename,states)
    character(*), intent(in) :: filename
    real(real64), allocatable, intent(out) :: states(:,:)
    integer, parameter :: capacity=100
    real(real64), parameter :: missing=-huge(1d0)
    real(real64) :: pressures(capacity),co2(capacity),h2o(capacity),co(capacity),soot(capacity)
    real(real64) :: tgas(capacity),tref(capacity)
    integer :: u,ios,n(7),i,j,k,l,m,t,r,row
    character(1024) :: message
    namelist /grid/ pressures,co2,h2o,co,soot,tgas,tref
    pressures=missing; co2=missing; h2o=missing; co=missing; soot=missing; tgas=missing; tref=missing
    open(newunit=u,file=filename,status='old',action='read',iostat=ios,iomsg=message)
    if (ios/=0) then
      write(error_unit,'(a)') trim(message)
      error stop 'Cannot open grid'
    endif
    read(u,nml=grid,iostat=ios,iomsg=message)
    close(u)
    if (ios/=0) then
      write(error_unit,'(a)') trim(message)
      error stop 'Cannot read grid namelist'
    endif
    n=[count(pressures/=missing),count(co2/=missing),count(h2o/=missing),count(co/=missing), &
       count(soot/=missing),count(tgas/=missing),count(tref/=missing)]
    if (any(n==0)) error stop 'Every grid axis must contain at least one value'
    if (product(int(n,int64))>1000000_int64) error stop 'Grid exceeds one million states'
    allocate(states(7,product(n)))
    row=0
    do i=1,n(1)
    do j=1,n(2)
    do k=1,n(3)
    do l=1,n(4)
    do m=1,n(5)
    do t=1,n(6)
    do r=1,n(7)
      row=row+1
      states(:,row)=[pressures(i),co2(j),h2o(k),co(l),soot(m),tgas(t),tref(r)]
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
  end subroutine
end module
