!!$**********************************************************************************
!!$                   SPECTRAL RADIATION CALCULATION SOFTWARE ~ SRCS ~
!!$**********************************************************************************
!!$
!!$                     COPYRIGHT (C) 2015 MICHAEL F. MODEST
!!$
!!$                           SCHOOL OF ENGINEERING                      
!!$                      UNIVERSITY OF CALIFORNIA AT MERCED
!!$
!!$	License: 
!!$		This file is a part of SRCS.
!!$		Used for non-commercial purposes only.
!!$		Not to be distributed without the permission of the copyright holder
!!$
!!$**********************************************************************************
module fskTableV4
!This module included the necessary subroutines for using NEW DEVELOPED FSK look-up table.
!It is firstly written by Chaojun Wang in Nov 2016.
implicit none

private

integer,parameter :: kFile=110
integer,parameter :: nPDB=34,nTgDB=28,nxco2DB=10,nxh2oDB=13,nxcoDB=6,nfvDB=4,NDB=32,Nq=32,Ndm=7
real(8),parameter :: PDB(nPDB)=(/0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.7d0,1.d0,2.d0,3.d0,4.d0,5.d0,6.d0,7.d0, &
                                 8.d0,9.d0,10.d0,11.d0,12.d0,13.d0,14.d0,15.d0,20.d0,25.d0,30.d0,35.d0, &
                                 40.d0,45.d0,50.d0,55.d0,60.d0,65.d0,70.d0,75.d0,80.d0/)
real(8),parameter :: xCO2DB(nxco2DB)=(/0.d0,0.01d0,0.02d0,0.03d0,0.04d0,0.05d0,0.25d0,0.5d0,0.75d0,1.d0/)
real(8),parameter :: xH2ODB(nxh2oDB)=(/0.d0,0.01d0,0.02d0,0.03d0,0.04d0,0.05d0,0.1d0,0.15d0,0.2d0,0.25d0,0.5d0,0.75d0,1.d0/)
real(8),parameter :: xCODB(nxcoDB)=(/0.d0,0.01d0,0.05d0,0.1d0,0.25d0,0.5d0/)
real(8),parameter :: fvDB(nfvDB)=(/0.d0,1.d-7,1.d-6,1d-5/)
real(8),parameter :: dT=100.d0,dTg=100.d0,cmin=1d-3,xmin=1d-7,fvmin=1d-9

TYPE kvalues
  REAL(4),dimension(Nq) :: values=0.d0
END TYPE kvalues

TYPE(kvalues),allocatable :: kq(:,:,:,:,:,:,:)
integer :: nP, nTg, nT, nfv, nxco2, nxh2o, nxco
integer :: iP,ifv,ixco2,ixh2o,ixco,iTg,iT,iq
real(8),dimension(:),allocatable :: P, Tg, T, fv, xCO2, xH2O, xCO
real(8) :: we,wP,wfv,wxco2,wxh2o,wxco,wTg,wT

  
public :: openFSKV4Database,closeFSKV4Database,fsk_v4_table,loadIntoMem,freeMem,get_ka,quadgen_
private :: coef_setup,maxID,detmLP,detmID
  
save

contains

integer function openFSKV4Database()

  real(8) :: Pm(2),Tmax,mole_frac_max(3),fvmax
  openFSKV4Database=0

  !Obtain k-values from V4 fsk look-up table
  Pm(1)=1.0
  Pm(2)=1.0 ! We need to parameterize this similar to LBL
  Tmax=3000
  mole_frac_max(1)=0.25 !CO2
  mole_frac_max(2)=0.25  !H2O
  mole_frac_max(3)=0.1 !CO
  fvmax=1.E-4
  
  call loadIntoMem(Pm,Tmax,mole_frac_max,fvmax)

  return
end function openFSKV4Database

integer function closeFSKV4Database()
  ! Nothing to done now.
  closeFSKV4Database=0
  return
end function closeFSKV4Database

subroutine fsk_v4_table (Mix_Info,Mix_Inforef,Nq, k, a, g, gFSK, wFSK, NoPDB)
use commonData, only : GasMixInfo
use fskTable, only: k_interp, a_interp

type(GasMixInfo),intent(in) :: Mix_Info,Mix_Inforef
integer, intent(in) :: Nq, NoPDB
real(8), dimension(Nq), intent(inout) :: k, a, g
real(8), dimension(NoPDB), intent(in) :: gFSK,wFSK
real(8), dimension(NoPDB) :: kpre, apre

call get_ka(Mix_Info,Mix_Inforef%T,gFSK,wFSK,kpre,apre)

call k_interp(NoPDB, gFSK, kpre,Nq, g, k)
call a_interp(NoPDB, gFSK, apre,Nq, g, a)

end subroutine fsk_v4_table

!--------------------------------------------------------------------------------------

subroutine loadIntoMem(Pm,Tmax,xmax,fvmax)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! function to open FSK look-up table
! currently table includes CO2,H2O and CO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use commonData, only : CO2, H2O, CO
use commonRoutines, only : strpad
use FSKDBPath, only : v4Path

implicit none

  real(8),intent(in) :: Pm(2),Tmax,xmax(3),fvmax

  ! local variables
  integer,parameter :: fK=101
  integer :: i,iPmax,iPmin,iTmax,id
  real(4) :: kDB(NDB+1)
  character(250) :: fNamek, fileNamek



  ! Determine index of minimum and maximum pressures
  call detmID(PDB,Pm(1),iPmin)
  call detmID(PDB,Pm(2),iPmax)
  if (iPmin==iPmax) then
    if ((Pm(2)-Pm(1))<=cmin) then
      if ((Pm(1)-PDB(iPmin))/(PDB(iPmin+1)-PDB(iPmin))<=cmin) then
        nP=1
      elseif ((Pm(1)-PDB(iPmin))/(PDB(iPmin+1)-PDB(iPmin))>=1.d0-cmin) then
        nP=1; iPmin=iPmin+1; iPmax=iPmin
      else
        nP=2; iPmax=iPmax+1
      endif
    else
      nP=2; iPmax=iPmax+1
    endif
  else
    iPmax=iPmax+1; nP=iPmax-iPmin+1
  endif
  if (allocated(P)) deallocate(P)
  allocate(P(nP))
  P(:)=PDB(iPmin:iPmax)

  ! Determine index of maximum temperature
  iTmax=min(nTgDB,int((Tmax-300.d0)/dTg)+2)
  nTg=iTmax; nT=nTg
  if (allocated(Tg)) deallocate(Tg)
  if (allocated(T)) deallocate(T)
  allocate(Tg(nTg), T(nT))
  do i=1,nTg
    Tg(i)=300.d0+(i-1)*dTg
    T(i)=300.d0+(i-1)*dT
  enddo

  ! Determine index of maximum mole fractions and soot volume fraction
  call maxID(xCO2DB,xmax(CO2),nxco2,.true.); nxco2=min(nxco2,nxco2DB)
  call maxID(xH2ODB,xmax(H2O),nxh2o,.true.); nxh2o=min(nxh2o,nxh2oDB)
  call maxID(xCODB,xmax(CO),nxco,.true.); nxco=min(nxco,nxcoDB)
  call maxID(fvDB,fvmax,nfv,.false.); nfv=min(nfv,nfvDB)
  if (allocated(xCO2)) deallocate(xCO2)
  if (allocated(xH2O)) deallocate(xH2O)
  if (allocated(xCO)) deallocate(xCO)
  if (allocated(fv)) deallocate(fv)
  allocate(xCO2(nxco2), xH2O(nxh2o), xCO(nxco), fv(nfv))
  xCO2(:)=xCO2DB(1:nxco2)
  xH2O(:)=xH2ODB(1:nxh2o)
  xCO(:)=xCODB(1:nxco)
  fv(:)=fvDB(1:nfv)

  ! Clear and allocate arrays
  if (allocated(kq)) deallocate(kq)
  allocate(kq(nP,nfv,nxco2,nxh2o,nxco,nTg,nT))

  ! Load k-values into arrays
  do iP=1,nP
    do ifv=1,nfv
      do ixco2=1,nxco2
        do ixh2o=1,nxh2o
          do ixco=1,nxco
            write(fNamek,91) P(iP),P(iP),xCO2(ixco2),xH2O(ixh2o),xCO(ixco),fv(ifv)
            call strpad(fNamek)
            fileNamek= ADJUSTL(ADJUSTR(v4Path)//ADJUSTL(fNamek))
            open(fK,file=fileNamek,form='unformatted',access='direct',recl=4*(NDB+1),status='old',action='read')
            do iTg=1,nTg
              id=nTgDB*(iTg-1)
              do iT=1,nT
                id=id+1
                read(fK,rec=id) (kDB(iq),iq=1,NDB+1)
                do iq=1,Nq
                  kq(iP,ifv,ixco2,ixh2o,ixco,iTg,iT)%values(iq)=kDB(1+(iq-1)*NDB/Nq)
                enddo
              enddo !iT
            enddo !iTg
            close(fK)
          enddo !ixco
        enddo !ixh2o
      enddo !ixco2
    enddo !ifv
  enddo !iP

  !write(*,*) 'Loading k-values >>> successful'

  91 FORMAT('/',F4.1,'/','P.',F4.1,'_CO2.',F4.2,'_H2O.',F4.2,'_CO.', F4.2,'_fv.',ES6.0,'_k.dat')

  return

end subroutine loadIntoMem

!-----------------------------------------------------------------------------------------------

subroutine freeMem()    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! function to close FSK look-up table
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
    
  integer :: ierr

  deallocate(P, Tg, T, fv, xCO2, xH2O, xCO, kq, stat=ierr)

  if (ierr==0) then
    write(*,*) 'Deallocating >>> successful'
  else
    write(*,*) 'Deallocating >>> failure'
  end if

  return

end subroutine freeMem

!------------------------------------------------------------------------------------------

subroutine get_ka(Mix_Info,Tref,gq,wq,k,a)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! subroutine to calculate the k-values and a-values for a given state
! by data interpolation using NEW DEVELOPED look-up table
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use commonData, only : GasMixInfo
use fskTable, only: k_interp
use libJCai, only: afun

implicit none
  
  TYPE(GasMixInfo),INTENT(IN) :: Mix_Info ! Local mixture information  
  real(8),intent(in) :: Tref ! reference temperature
  real(8),dimension(Nq),intent(in) :: gq,wq !quadrature points
  real(8),dimension(Nq),intent(out) :: k,a ! k-values and a-values

  real(8),dimension(Nq) :: kll,gql ! Note that gql is useless

  ! get k-values
  call get_k(Mix_Info,Tref,k)

  ! get a-values
  call get_k(Mix_Info,Mix_Info%T,kll)
  call afun(gq, k, gq, kll, gql, wq, a)

  return

end subroutine get_ka

!------------------------------------------------------------------------------------------

subroutine get_k(Mix_Info,Tref,k)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! subroutine to calculate the k-values for a given state
! by data interpolation using NEW DEVELOPED look-up table
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use commonData, only : GasMixInfo

implicit none
  
  TYPE(GasMixInfo),INTENT(IN) :: Mix_Info ! Local mixture information  
  real(8),intent(in) :: Tref ! reference temperature
  real(8),dimension(Nq),intent(out) :: k ! k-values

  ! local variables
  integer :: grid(Ndm,2),lp(Ndm)
  real(8) :: coef(Ndm)
    
  ! initialization of output data
  k=0.d0; grid=0; coef=0.d0

  call coef_setup(Mix_Info,Tref,grid,coef)
  call detmLP(coef,lp)

  do iP=1,lp(1)
    if (iP==1) then; wP=1.d0-coef(1)
    else; wP=coef(1); endif
    do ifv=1,lp(2)
      if (ifv==1) then; wfv=1.d0-coef(2)
      else; wfv=coef(2); endif
      do ixco2=1,lp(3)
        if (ixco2==1) then; wxco2=1.d0-coef(3)
        else; wxco2=coef(3); endif
        do ixh2o=1,lp(4)
          if (ixh2o==1) then; wxh2o=1.d0-coef(4)
          else; wxh2o=coef(4); endif
          do ixco=1,lp(5)
            if (ixco==1) then; wxco=1.d0-coef(5)
            else; wxco=coef(5); endif
            do iTg=1,lp(6)
              if (iTg==1) then; wTg=1.d0-coef(6)
              else; wTg=coef(6); endif
              do iT=1,lp(7)
                if (iT==1) then; wT=1.d0-coef(7)
                else; wT=coef(7); endif

                !calculate weight
                we=wP*wfv*wxco2*wxh2o*wxco*wTg*wT

                !calculate k-values
                k(1:Nq)=k(1:Nq)+we*kq(grid(1,iP),grid(2,ifv),grid(3,ixco2),grid(4,ixh2o),grid(5,ixco),grid(6,iTg),grid(7,iT)) &
                                   %values(1:Nq)

              enddo !iT
            enddo !iTg
          enddo !ixco
        enddo !ixh2o
      enddo !ixco2
    enddo !ifv
  enddo !iP

  return

end subroutine get_k

!-----------------------------------------------------------------------------------------

subroutine coef_setup(Mix_Info,Tref,grid,coef)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! find the interpolation grid and coefficients in composition vaiables space
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use commonData, only : GasMixInfo

IMPLICIT NONE

  TYPE(GasMixInfo),INTENT(IN):: Mix_Info ! gas mixture information
  REAL(8),INTENT(IN) ::  Tref ! reference temperature
  INTEGER,INTENT(OUT) :: grid(Ndm,2) ! grid
  REAL(8),INTENT(OUT) :: coef(Ndm) ! coefficients

  grid=0; coef=0.d0

  ! Pressure index
  if (nP==1) then
    grid(1,1)=1; coef(1)=0.d0
  else
    call detmID(P,Mix_Info%P,grid(1,1))
    grid(1,2)=grid(1,1)+1
    coef(1)=(Mix_Info%P-P(grid(1,1)))/(P(grid(1,2))-P(grid(1,1)))
  endif

  ! fv index
  if (nfv==1) then
    grid(2,1)=1; coef(2)=0.d0
  else
    call detmID(fv,Mix_Info%fv,grid(2,1))
    grid(2,2)=grid(2,1)+1
    coef(2)=(Mix_Info%fv-fv(grid(2,1)))/(fv(grid(2,2))-fv(grid(2,1)))
  endif

  ! CO2 index
  if (nxco2==1) then
    grid(3,1)=1; coef(3)=0.d0
  else
    call detmID(xCO2,Mix_Info%xCO2,grid(3,1))
    grid(3,2)= grid(3,1)+1
    coef(3)=(Mix_Info%xCO2-xCO2(grid(3,1)))/(xCO2(grid(3,2))-xCO2(grid(3,1)))
  endif

  ! H2O index
  if (nxh2o==1) then
    grid(4,1)=1; coef(4)=0.d0
  else
    call detmID(xH2O,Mix_Info%xH2O,grid(4,1))
    grid(4,2)= grid(4,1)+1
    coef(4)=(Mix_Info%xH2O-xH2O(grid(4,1)))/(xH2O(grid(4,2))-xH2O(grid(4,1)))
  endif

  ! CO index
  if (nxco==1) then
    grid(5,1)=1; coef(5)=0.d0
  else
    call detmID(xCO,Mix_Info%xCO,grid(5,1))
    grid(5,2)= grid(5,1)+1
    coef(5)=(Mix_Info%xCO-xCO(grid(5,1)))/(xCO(grid(5,2))-xCO(grid(5,1)))
  endif

  ! Local temperature index
  grid(6,1)= 1+INT((Mix_Info%T-Tg(1))/dTg)
  grid(6,1)= MAX(grid(6,1),1); grid(6,1)= MIN(grid(6,1),nTg-1)
  grid(6,2)= grid(6,1)+1
  coef(6)=(Mix_Info%T-Tg(grid(6,1)))/dTg

  ! Reference temperature index
  grid(7,1)= 1+INT((Tref-T(1))/dT)
  grid(7,1)= MAX(grid(7,1),1); grid(7,1)= MIN(grid(7,1),nT-1)
  grid(7,2)= grid(7,1)+1
  coef(7)=(Tref-T(grid(7,1)))/dT

  return

end subroutine coef_setup

!------------------------------------------------------------------------------------------

subroutine maxID(arr,val,id,isGas)

implicit none

  real(8),dimension(:),intent(in) :: arr
  real(8),intent(in) :: val
  integer,intent(out) :: id
  logical,intent(in) :: isGas

  integer :: i,imax
  real(8) :: valmin

  if (isGas) then
    valmin=xmin
  else
    valmin=fvmin
  endif
    
  if (val<=valmin) then
    id=1
  else
    call detmID(arr,val,id); id=id+1
  endif

  return

end subroutine maxID

!-----------------------------------------------------------------------------------------

subroutine detmLP(arr,val)
  implicit none

  real(8),dimension(:),intent(in) :: arr
  integer,dimension(size(arr)),intent(out) :: val

  integer :: i

  do i=1,size(arr)
    if (arr(i)<=cmin) then; val(i)=1
    else; val(i)=2; endif
  enddo

end subroutine detmLP

!------------------------------------------------------------------------------------------

subroutine detmID(arr,val,id)
  implicit none

  real(8),dimension(:),intent(in) :: arr
  real(8),intent(in) :: val
  integer,intent(out) :: id

  integer :: lef,rig,mid

  if (val>=arr(size(arr))) then
    id=size(arr)-1
  elseif (val<=arr(1)) then
    id=1
  else
    lef=1; rig=size(arr); mid=(lef+rig)/2
    do while (rig>lef+1)
      IF (val< arr(mid)) THEN
        rig= mid
      ELSE
        lef= mid
      END IF
      mid= (lef+rig)/2
    END DO
    id=lef
  endif

end subroutine detmID

!------------------------------------------------------------------------------------------

SUBROUTINE quadgen_(Cheb2,g,w,NoP,alpha) !quadrature closed at 0
  use quadlib, only : gaucheb2
  IMPLICIT NONE
  LOGICAL,INTENT(IN)    :: Cheb2
  INTEGER,INTENT(IN)    :: NoP
  REAL(8),INTENT(IN),optional    :: alpha  ! transformation factor
  REAL(8),INTENT(OUT)   :: g(NoP),w(NoP)
  !local variables
  REAL(8) :: gg(2*NoP-1),ww(2*NoP-1),dummy
  INTEGER(2) :: i

    CALL gaucheb2(gg,ww,2*NoP-1)
    g= -gg(NoP:2*NoP-1)
    w= ww(NoP:2*NoP-1)
    g(1)=abs(g(1)); w(1)=w(1)/2.d0
    IF ((.not.Cheb2).and.(.not. PRESENT(alpha))) THEN                ! transformation on Gauss-Chebyshev2 quadrature
      DO i=1, NoP
        w(i)= 1.5d0*SQRT(1.d0-g(i))*w(i)
        g(i)= 1.d0-(1.d0-g(i))**1.5
      ENDDO
    ELSEIF((.not.Cheb2).and.(PRESENT(alpha))) THEN
      DO i= 1, NoP
        w(i)= alpha*w(i)*(1.d0-g(i))**(alpha-1.d0)
        g(i)= 1.d0-(1.d0-g(i))**alpha
      ENDDO
    ENDIF
    dummy= SUM(w(1:NoP))
    w(1:NoP)=w(1:NoP)/dummy

END SUBROUTINE quadgen_

!-----------------------------------------------------------------------------

end module fskTableV4
