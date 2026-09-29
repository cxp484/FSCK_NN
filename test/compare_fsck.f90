! Standalone comparison: Fortran calls to ONNX Runtime and original V4 interpolation.
program compare_fsck
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic
  use comparison_grid, only: read_grid
  use nn_model
  use commonData, only: GasMixInfo
  use FSKDBPath, only: v4path
  use fskTableV4, only: loadIntoMem, get_ka, freeMem
  use libJCai, only: afun
  use quadlib, only: quadgen2
  implicit none
  character(1024) :: model_file,states_file,output_file,dbpath
  real(real64), allocatable :: states(:,:)
  real(real64) :: g(32),w(32),kn(32),kl(32),kd(32),an(32),ad(32),local(7)
  real(real64) :: pm(2),xmax(3),tmax,fvmax
  type(GasMixInfo) :: mix
  integer :: u,v,n,i,q,mono
  if (command_argument_count()/=4) error stop 'Usage: compare_fsck model.onnx grid.nml database output.csv'
  call get_command_argument(1,model_file)
  call get_command_argument(2,states_file)
  call get_command_argument(3,dbpath)
  call get_command_argument(4,output_file)
  call read_model(trim(model_file))
  call read_grid(trim(states_file),states)
  n=size(states,2)
  if (.not.all(ieee_is_finite(states))) error stop 'Nonfinite state'
  if (any(states(1,:)<.1d0).or.any(states(1,:)>80d0)) error stop 'Pressure out of database range'
  if (any(states(2:4,:)<0d0).or.any(states(2:3,:)>1d0).or.any(states(4,:)>.5d0)) error stop 'Gas out of range'
  if (any(sum(states(2:4,:),dim=1)>1d0+1d-12)) error stop 'Mole fractions exceed unity'
  if (any(states(5,:)<0d0).or.any(states(5,:)>1d-5)) error stop 'Soot out of range'
  if (any(states(6:7,:)<300d0).or.any(states(6:7,:)>3000d0)) error stop 'Temperature out of range'
  ! Use the unmodified loader, including its original bracketing and tolerances.
  if (len_trim(dbpath)>150) error stop 'Database path exceeds original loader filename capacity'
  v4path=trim(dbpath)
  pm=[minval(states(1,:)),maxval(states(1,:))]
  xmax=maxval(states(2:4,:),dim=2)
  tmax=maxval(states(6:7,:)); fvmax=maxval(states(5,:))
  call loadIntoMem(pm,tmax,xmax,fvmax)
  call quadgen2(.false.,g,w,32,2d0)
  open(newunit=v,file=trim(output_file),status='replace',action='write')
  write(v,'(a)') 'state_id,q,P_atm,xCO2,xH2O,xCO,fv,Tgas_K,Tref_K,g,w,k_db,k_nn,a_db,a_nn,nn_monotone'
  do i=1,n
    mix%P=states(1,i); mix%xCO2=states(2,i); mix%xH2O=states(3,i)
    mix%xCO=states(4,i); mix%fv=states(5,i); mix%T=states(6,i)
    mix%xCH4=0d0; mix%xC2H4=0d0
    call get_ka(mix,states(7,i),g,w,kd,ad)
    call predict(states(:,i),kn)
    local=states(:,i); local(7)=local(6)
    call predict(local,kl)
    if (.not.all(ieee_is_finite(kn)).or..not.all(ieee_is_finite(kl))) error stop 'Nonfinite NN kappa'
    mono=0
    if (all(kn(2:)>=kn(:31)).and.all(kl(2:)>=kl(:31))) mono=1
    ! Exactly the same two-distribution construction as get_ka.
    ! Nonmonotone curves are retained and flagged: a is diagnostic in those cases.
    call afun(g,kn,g,kl,g,w,an)
    do q=1,32
      write(v,'(i0,",",i0,13(",",es24.16e3),",",i0)') &
        i,q,states(:,i),g(q),w(q),kd(q),kn(q),ad(q),an(q),mono
    enddo
  enddo
  close(v)
  call freeMem()
  call close_model()
  open(newunit=u,file=trim(output_file)//'.meta',status='replace',action='write')
  write(u,'(a)') 'model='//trim(model_file), 'grid='//trim(states_file), 'database='//trim(dbpath)
  close(u)
  write(*,'(a,i0,a)') 'Compared ',n,' states.'
end program
