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
module fskTable
!This module included the necessary subroutines for using FSK look-up table.
!It is firstly written by Wenjun Ge in Oct 2014 and modified by Chaojun Wang
!in Apr 2015.
implicit none

private

TYPE kavalues
  REAL(4),dimension(:),pointer :: values => null()
END TYPE kavalues

integer,parameter :: kFile=110,aFile=120
integer,parameter :: nP=34,nT=28,nTp=28,nxco2=10,nxh2o=13,nxco=6,nq=32,num=2
integer,parameter :: xincco2=39200,xinch2o=3920,xincco=784,Tinc=28
integer,parameter :: idco2(0:nxco2-1)=(/0,77,154,231,308,385,462,539,610,675/)
integer,parameter :: idh2o(4,0:nxh2o-1)=reshape((/0,0,0,0,6,6,6,5,12,12,12,10,18,18,18,15, &
						  24,24,24,20,30,30,30,25,36,36,36,30, &
						  42,42,42,35,48,48,48,40,54,54,54,44, &
						  60,60,60,47,66,66,65,47,72,71,65,47/),(/4,nxh2o/))
TYPE(kavalues),dimension(:,:,:,:,:,:),allocatable :: kq
TYPE(kavalues),dimension(:,:,:,:,:,:),allocatable :: aq
real(8),parameter :: xCO2(nxco2)=(/0.d0,0.01d0,0.02d0,0.03d0,0.04d0,0.05d0,0.25d0,0.5d0,0.75d0,1.d0/) 
real(8),parameter :: xH2O(nxh2o)=(/0.d0,0.01d0,0.02d0,0.03d0,0.04d0,0.05d0,0.1d0,0.15d0,0.2d0,0.25d0,0.5d0,0.75d0,1.d0/)
real(8),parameter :: xCO(nxco)=(/0.d0,0.01d0,0.05d0,0.1d0,0.25d0,0.5d0/)
real(8),parameter :: P(nP)=(/0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.7d0,1.d0,2.d0,3.d0,4.d0,5.d0,6.d0,7.d0, &
			     8.d0,9.d0,10.d0,11.d0,12.d0,13.d0,14.d0,15.d0,20.d0,25.d0,30.d0,35.d0, &
			     40.d0,45.d0,50.d0,55.d0,60.d0,65.d0,70.d0,75.d0,80.d0/)
real(8),parameter :: T(nT)= (/300.d0,400.d0,500.d0,600.d0,700.d0,800.d0,900.d0,1000.d0,1100.d0, &
     1200.d0,1300.d0,1400.d0,1500.d0,1600.d0,1700.d0,1800.d0,1900.d0,2000.d0,2100.d0,2200.d0, &
     2300.d0,2400.d0,2500.d0,2600.d0,2700.d0,2800.d0,2900.d0,3000.d0/)
real(8),parameter :: Tp(nTp)= (/300.d0,400.d0,500.d0,600.d0,700.d0,800.d0,900.d0,1000.d0,1100.d0, &
     1200.d0,1300.d0,1400.d0,1500.d0,1600.d0,1700.d0,1800.d0,1900.d0,2000.d0,2100.d0,2200.d0, &
     2300.d0,2400.d0,2500.d0,2600.d0,2700.d0,2800.d0,2900.d0,3000.d0/)
real(8),parameter :: dxco2_1=0.01d0,dxco2_2=0.2d0,dxco2_3=0.25d0,dxh2o_1=0.01d0,dxh2o_2=0.05d0,dxh2o_3=0.25d0, &
		     dxco_1=0.01d0,dxco_2=0.04d0,dxco_3=0.05d0,dxco_4=0.15d0,dxco_5=0.25d0
real(8),parameter :: dT=100.d0,dTp=100.d0
real(8) :: weight,weP,weTp,weT,wexco2,wexh2o,wexco
real(8) :: wP,wTp,wT,wxco2,wxh2o,wxco
  
public :: openDatabase,closeDatabase,kgaInterp,kgaInterpM,k_interp,a_interp,fsk_table
private :: gridPos,get_k,get_a,get_k_m,get_a_m
  
save

contains

!--------------------------------------------------------------------------------------

integer function openDatabase(tableNum)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! function to open FSK look-up table
! currently table includes CO2,H2O and CO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use FSKDBPath, only : ktPath

implicit none

  integer, intent(in) :: tableNum !table number
  ! local variables
  integer :: ierr1,ierr2
  character(512) :: DBk,DBa,fk,fa

  write(*,*)"ktpath=",ktPath

  openDatabase=0
  ! open data files
  if (tableNum == 1) then
    allocate(kq(nP,nTp,nT,nxco2-3,nxh2o-3,nxco-1),aq(nP,nTp,nT,nxco2-3,nxh2o-3,nxco-1))
    fk='FSK_kTable_Part.dat'; fa='FSK_aTable_Part.dat'
    DBk=ADJUSTL(ADJUSTR(ktPath)//ADJUSTL(fk))
    DBa=ADJUSTL(ADJUSTR(ktPath)//ADJUSTL(fa))
    open(unit=kFile, file=DBk, form = 'unformatted', access='direct', &
         recl=Nq*4, status='old', action = 'read', iostat=ierr1)
    open(unit=aFile, file=DBa, form = 'unformatted', access='direct', &
         recl=Nq*4, status='old', action = 'read', iostat=ierr2)
    if (ierr1==0 .and. ierr2==0) then
      write(*,*) 'Opening FSK Look-up Part Table >>> successful'
      openDatabase=0
    else
      write(*,*) 'Opening FSK Look-up Part Table >>> failure'
      openDatabase=1
    end if
  elseif (tableNum == 2) then
    allocate(kq(nP,nTp,nT,nxco2,nxh2o,nxco),aq(nP,nTp,nT,nxco2,nxh2o,nxco))
    fk='FSK_kTable.dat'; fa='FSK_aTable.dat'
    DBk=ADJUSTL(ADJUSTR(ktPath)//ADJUSTL(fk))
    DBa=ADJUSTL(ADJUSTR(ktPath)//ADJUSTL(fa))
    open(unit=kFile, file=DBk, form = 'unformatted', access='direct', &
         recl=Nq*4, status='old', action = 'read', iostat=ierr1)
    open(unit=aFile, file=DBa, form = 'unformatted', access='direct', &
         recl=Nq*4, status='old', action = 'read', iostat=ierr2)
    if (ierr1==0 .and. ierr2==0) then
      write(*,*) 'Opening FSK Look-up Table >>> successful'
      openDatabase=0
    else
      write(*,*) 'Opening FSK Look-up Table >>> failure'
      openDatabase=1
    end if

  endif

  return

end function openDatabase

!-----------------------------------------------------------------------------------------------

integer function closeDatabase()    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! function to close FSK look-up table
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
implicit none
    
  integer :: ierr

  close(kFile,iostat=ierr)
  if (ierr==0) then
    write(*,*) 'closing k-Table >>> successful'
    closeDatabase=0
  else
    write(*,*) 'closing k-Table >>> failure'
    closeDatabase=1
  end if
  close(aFile,iostat=ierr)
  if (ierr==0) then
    write(*,*) 'closing a-Table >>> successful'
    closeDatabase=0
  else
    write(*,*) 'closing a-Table >>> failure'
    closeDatabase=1
  end if

  return

end function closeDatabase

!------------------------------------------------------------------------------------------

subroutine kgaInterp(Mix_Info,Mix_Inforef,Tref,k,a,gq)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! subroutine to calculate the k-values and a-values for a given state
! by data interpolation using XQtrTable
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use commonData, only : GasMixInfo

implicit none
  
  TYPE(GasMixInfo),INTENT(IN) :: Mix_Info,Mix_Inforef !Local and reference mixture information  
  real(8),intent(in) :: Tref !reference temperature
  real(8),dimension(Nq),intent(in) :: gq !quadrature points
  real(8),dimension(Nq),intent(out) :: k,a !k-values and a-values
  ! local variables
  integer :: i,iq,grid(6,num)
  real(8),dimension(Nq) :: kll,krr,krl,arl,gql
  ! kll: k with absc eval at local thermo state, weighted by Planck fun at local temperature
  ! krl: k with absc eval at reference thermo state, weighted by Planck fun at local temperature
  ! krr: k with absc eval at ref thermo state, weighted by Planck fun at ref temperature 
  real(8) :: Val=10.0d0**(-6.0d0)   
    
  ! initialization of output data
  k=0.d0; a=1.d0; grid=0; arl=1.d0; kll=0.d0; krr=0.d0; krl=0.d0

  call gridPos(Mix_Info,Mix_Info%T,grid) 
  call get_k(grid,Mix_Info,Mix_Info%T,kll)

  call gridPos(Mix_Inforef,Mix_Inforef%T,grid)
  call get_k(grid,Mix_Inforef,Mix_Inforef%T,krr)

  call gridPos(Mix_Inforef,Mix_Info%T,grid)
  call get_k(grid,Mix_Inforef,Mix_Info%T,krl)

  call k_interp(Nq,krl,gq,Nq,krr,gql)
  call k_interp(Nq,gq,kll,Nq,gql,k)

  !a function
  call get_a(grid,Mix_Inforef,Mix_Info%T,arl)
  call a_interp(Nq,gq,arl,Nq,gql,a)

end subroutine kgaInterp

!-----------------------------------------------------------------------------------------

subroutine kgaInterpM(Mix_Info,Mix_Inforef,Tref,k,a,gq)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! subroutine to calculate the k-values and a-values for a given state
! by data interpolation using FullTable
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use commonData, only : GasMixInfo

implicit none
  
  TYPE(GasMixInfo),INTENT(IN) :: Mix_Info,Mix_Inforef  
  real(8),intent(in) :: Tref
  real(8),dimension(Nq),intent(in) :: gq
  real(8),dimension(Nq),intent(out) :: k,a
 
  integer :: i,iq,grid(6,num)
  real(8),dimension(Nq) :: kll,krr,krl,arl,gql
  real(8) :: Val=10.0d0**(-6.0d0)   
    
  ! initialization of output data
  k=0.d0; a=1.d0; grid=0; arl=1.d0; kll=0.d0

  call gridPos(Mix_Info,Mix_Info%T,grid) 
  call get_k_m(grid,Mix_Info,Mix_Info%T,kll)

  call gridPos(Mix_Inforef,Mix_Inforef%T,grid)
  call get_k_m(grid,Mix_Inforef,Mix_Inforef%T,krr)

  call gridPos(Mix_Inforef,Mix_Info%T,grid)
  call get_k_m(grid,Mix_Inforef,Mix_Info%T,krl)

  call k_interp(Nq,krl,gq,Nq,krr,gql)
  call k_interp(Nq,gq,kll,Nq,gql,k)

  !a function
  call get_a_m(grid,Mix_Inforef,Mix_Info%T,arl)
  call a_interp(Nq,gq,arl,Nq,gql,a)

end subroutine kgaInterpM

!-----------------------------------------------------------------------------------------

subroutine gridPos(Mix_Info,Trad,grid)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! find the interpolation grid in composition vaiables space
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use commonData, only : GasMixInfo

IMPLICIT NONE

  TYPE(GasMixInfo),INTENT(IN):: Mix_Info !gas mixture information
  real(8),intent(in) ::  Trad !Planck temperature
  INTEGER,INTENT(OUT) :: grid(6,num) !grid
  ! local variables
  INTEGER(2):: iPa(num),iTa(num),iTpa(num),ico2a(num),ih2oa(num),icoa(num)
  INTEGER(2):: iP,i,j 

  !Locate pressure grid
  if (Mix_Info%P < P(6)) then !0.7
    i=int((Mix_Info%P-P(1))/0.1d0)+1
  elseif (Mix_Info%P < P(7)) then !1.0
    i=6
  elseif (Mix_Info%P < P(21)) then !15.0
    i=int(Mix_Info%P-P(7))+7
  elseif (Mix_Info%P <= P(34)) then !80.0
    i=int((Mix_Info%P-P(21))/5.d0)+21
  endif
  if (i>nP-1) i=nP-1
  iPa(1)=i; iPa(2)=i+1

  ! Locate Planck function temperature grid
  i=int((Trad-Tp(1))/dTp)+1
  IF (i<1) THEN;i=1
  ELSEIF(i>nTp-1) THEN;i=nTp-1
  ENDIF
  iTpa(1)= i; iTpa(2)= i+1

  ! Locate composition temperature grid
  i=int((Mix_Info%T-T(1))/dT)+1
  IF (i<1) THEN;i=1
  ELSEIF(i>nT-1) THEN;i=nT-1
  ENDIF
  iTa(1)= i; iTa(2)= i+1

  ! Locate co2 mole fraction grid
  if (Mix_Info%xCO2 < xCO2(6)) then
    i=int(Mix_Info%xCO2/dxco2_1)+1
  elseif (Mix_Info%xCO2 < xCO2(nxco2)) then
    i=int(Mix_Info%xCO2/dxco2_3)+6
  elseif (Mix_Info%xCO2 == xCO2(nxco2)) then
    i=9
  endif
  ico2a(1)= i; ico2a(2)= i+1

  ! Locate h2o mole fraction grid
  if (Mix_Info%xH2O < xH2O(6)) then
    i=int(Mix_Info%xH2O/dxh2o_1)+1
  elseif (Mix_Info%xH2O < xH2O(10)) then
    i=int(Mix_Info%xH2O/dxh2o_2)+5
  elseif (Mix_Info%xH2O < xH2O(nxh2o)) then
    i=int(Mix_Info%xH2O/dxh2o_3)+9
  elseif (Mix_Info%xH2O == xH2O(nxh2o)) then
    i=12
  endif
  ih2oa(1)= i; ih2oa(2)= i+1

  ! Locate co mole fraction grid
  if (Mix_Info%xCO < xCO(2)) then
    i=1
  elseif (Mix_Info%xCO < xCO(3)) then
    i=2
  elseif (Mix_Info%xCO < xCO(4)) then
    i=3
  elseif (Mix_Info%xCO < xCO(5)) then
    i=4
  else
    i=5
  endif
  icoa(1)= i; icoa(2)= i+1

  grid(1,:)= iPa; grid(2,:)= iTpa; grid(3,:)= iTa
  grid(4,:)= ico2a; grid(5,:)= ih2oa; grid(6,:)= icoa  

end subroutine gridPos

!-----------------------------------------------------------------------------------------

subroutine get_k(grid,Mix_Info,Trad,final)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Get k-values directly from XQtrTable and then carry on linear interpolation 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use commonData, only : GasMixInfo

  IMPLICIT NONE
  TYPE(GasMixInfo),INTENT(IN):: Mix_Info !gas mixture information
  INTEGER,INTENT(IN) :: grid(6,num) !grid
  real(8),intent(in) ::  Trad !Planck temperature
  real(8),dimension(Nq),intent(out) :: final !output k-values
  ! local variables
  integer,parameter :: Pinc=273616
  integer :: i,j,k,l,m,n,id,iq

  final=0.d0
  wP=(Mix_Info%P-P(grid(1,1)))/(P(grid(1,2))-P(grid(1,1)))
  wTp=(Trad-Tp(grid(2,1)))/dTp
  wT=(Mix_Info%T-T(grid(3,1)))/dT
  if (Mix_Info%xCO2 < xCO2(6)) then
    wxco2=(Mix_Info%xCO2-xCO2(grid(4,1)))/dxco2_1
  else
    wxco2=(Mix_Info%xCO2-xCO2(grid(4,1)))/dxco2_2
  endif
  if (Mix_Info%xH2O < xH2O(6)) then
    wxh2o=(Mix_Info%xH2O-xH2O(grid(5,1)))/dxh2o_1
  else
    wxh2o=(Mix_Info%xH2O-xH2O(grid(5,1)))/dxh2o_2
  endif
  if (Mix_Info%xCO < xCO(2)) then
    wxco=(Mix_Info%xCO-xCO(grid(6,1)))/dxco_1
  elseif (Mix_Info%xCO < xCO(3)) then
    wxco=(Mix_Info%xCO-xCO(grid(6,1)))/dxco_2
  elseif (Mix_Info%xCO < xCO(4)) then
    wxco=(Mix_Info%xCO-xCO(grid(6,1)))/dxco_3
  else
    wxco=(Mix_Info%xCO-xCO(grid(6,1)))/dxco_4
  endif

  do i=1,num
    if (i==1) then; weP=1.d0-wP;
    else; weP=wP
    endif
    do j=1,num
      if (j==1) then; weTp=1.d0-wTp;
      else; weTp=wTp
      endif
      do k=1,num
        if (k==1) then; weT=1.d0-wT;
        else; weT=wT
        endif
	do l=1,num
          if (l==1) then; wexco2=1.d0-wxco2;
          else; wexco2=wxco2
          endif
	  do m=1,num
            if (m==1) then; wexh2o=1.d0-wxh2o;
            else; wexh2o=wxh2o
            endif
	    do n=1,num
              if (n==1) then; wexco=1.d0-wxco;
              else; wexco=wxco
              endif
	      !calculate weight
	      weight=weP*weTp*weT*wexco2*wexh2o*wexco
	      !load k-values into memory
	      if (associated(kq(grid(1,i),grid(2,j),grid(3,k), &
				grid(4,l),grid(5,m),grid(6,n))%values)) then
	        continue
	      else
	        allocate(kq(grid(1,i),grid(2,j),grid(3,k),grid(4,l),grid(5,m),grid(6,n))%values(Nq))
		if (grid(4,l)==1 .and. grid(5,m)==1 .and. grid(6,n)==1) then
		  kq(grid(1,i),grid(2,j),grid(3,k),grid(4,l),grid(5,m),grid(6,n))%values=0.d0
		else
		  if (grid(4,l)==1) then
		    if (grid(5,m)==1) then
		      id=Pinc*(grid(1,i)-1)+xincco*(grid(6,n)-2)+Tinc*(grid(3,k)-1)+grid(2,j)
		    else
		      id=Pinc*(grid(1,i)-1)+xinch2o*(grid(5,m)-1)-xincco+ &
		         xincco*(grid(6,n)-1)+Tinc*(grid(3,k)-1)+grid(2,j)
		    endif
	          elseif (grid(4,l)==2) then
		    id=Pinc*(grid(1,i)-1)+(xincco2-xincco)*(grid(4,l)-1)+xinch2o*(grid(5,m)-1)+ &
		       xincco*(grid(6,n)-1)+Tinc*(grid(3,k)-1)+grid(2,j)
	          else !elseif (grid(4,l)==3) then
		    id=Pinc*(grid(1,i)-1)+xincco2*(grid(4,l)-1)-xincco+xinch2o*(grid(5,m)-1)+ &
		       xincco*(grid(6,n)-1)+Tinc*(grid(3,k)-1)+grid(2,j)
	          endif
	          read(kFile,rec=id) &
		    (kq(grid(1,i),grid(2,j),grid(3,k),grid(4,l),grid(5,m),grid(6,n))%values(iq),iq=1,Nq)
		endif
	      endif
	      !calculate k-distributions
	      final(1:Nq)=final(1:Nq)+weight*kq(grid(1,i),grid(2,j),grid(3,k), &
						grid(4,l),grid(5,m),grid(6,n))%values(1:Nq)
	    enddo !n
	  enddo !m
	enddo !l
      enddo !k
    enddo !j
  enddo !i

  return

end subroutine get_k

!-----------------------------------------------------------------------------------------

subroutine get_k_m(grid,Mix_Info,Trad,final)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Get k-values directly from FullTable and then carry on linear interpolation 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use commonData, only : GasMixInfo

  IMPLICIT NONE
  TYPE(GasMixInfo),INTENT(IN):: Mix_Info
  INTEGER,INTENT(IN) :: grid(6,num)
  real(8),intent(in) ::  Trad
  real(8),dimension(Nq),intent(out) :: final
  integer,parameter :: Pinc=565264
  integer :: i,j,k,l,m,n,id,iq

  final=0.d0
  wP=(Mix_Info%P-P(grid(1,1)))/(P(grid(1,2))-P(grid(1,1)))
  wTp=(Trad-Tp(grid(2,1)))/dTp
  wT=(Mix_Info%T-T(grid(3,1)))/dT
  if (Mix_Info%xCO2 < xCO2(6)) then
    wxco2=(Mix_Info%xCO2-xCO2(grid(4,1)))/dxco2_1
  elseif (Mix_Info%xCO2 < xCO2(7)) then
    wxco2=(Mix_Info%xCO2-xCO2(grid(4,1)))/dxco2_2
  elseif (Mix_Info%xCO2 <= xCO2(nxco2)) then
    wxco2=(Mix_Info%xCO2-xCO2(grid(4,1)))/dxco2_3
  endif
  if (Mix_Info%xH2O < xH2O(6)) then
    wxh2o=(Mix_Info%xH2O-xH2O(grid(5,1)))/dxh2o_1
  elseif (Mix_Info%xH2O < xH2O(10)) then
    wxh2o=(Mix_Info%xH2O-xH2O(grid(5,1)))/dxh2o_2
  elseif (Mix_Info%xH2O <= xH2O(nxh2o)) then
    wxh2o=(Mix_Info%xH2O-xH2O(grid(5,1)))/dxh2o_3
  endif
  if (Mix_Info%xCO < xCO(2)) then
    wxco=(Mix_Info%xCO-xCO(grid(6,1)))/dxco_1
  elseif (Mix_Info%xCO < xCO(3)) then
    wxco=(Mix_Info%xCO-xCO(grid(6,1)))/dxco_2
  elseif (Mix_Info%xCO < xCO(4)) then
    wxco=(Mix_Info%xCO-xCO(grid(6,1)))/dxco_3
  elseif (Mix_Info%xCO < xCO(5)) then
    wxco=(Mix_Info%xCO-xCO(grid(6,1)))/dxco_4
  else
    wxco=(Mix_Info%xCO-xCO(grid(6,1)))/dxco_5
  endif

  do i=1,num
    if (i==1) then; weP=1.d0-wP;
    else; weP=wP
    endif
    do j=1,num
      if (j==1) then; weTp=1.d0-wTp;
      else; weTp=wTp
      endif
      do k=1,num
        if (k==1) then; weT=1.d0-wT;
        else; weT=wT
        endif
	do l=1,num
          if (l==1) then; wexco2=1.d0-wxco2;
          else; wexco2=wxco2
          endif
	  do m=1,num
            if (m==1) then; wexh2o=1.d0-wxh2o;
            else; wexh2o=wxh2o
            endif
	    do n=1,num
              if (n==1) then; wexco=1.d0-wxco;
              else; wexco=wxco
              endif
	      !calculate weight
	      weight=weP*weTp*weT*wexco2*wexh2o*wexco
	      !load k-values into memory
	      if (associated(kq(grid(1,i),grid(2,j),grid(3,k), &
				grid(4,l),grid(5,m),grid(6,n))%values)) then
	        continue
	      else
	        allocate(kq(grid(1,i),grid(2,j),grid(3,k),grid(4,l),grid(5,m),grid(6,n))%values(Nq))
		if (grid(4,l)==1 .and. grid(5,m)==1 .and. grid(6,n)==1) then
		  kq(grid(1,i),grid(2,j),grid(3,k),grid(4,l),grid(5,m),grid(6,n))%values=0.d0
		else
		  if (grid(4,l)<8) then
		    if (grid(4,l)==1) then
		      if (grid(5,m)==1) then
			id=Pinc*(grid(1,i)-1)+xincco*(grid(6,n)-2)+Tinc*(grid(3,k)-1)+grid(2,j)
		      else
		        id=Pinc*(grid(1,i)-1)+xincco*(idh2o(1,(grid(5,m)-1))-1)+ &
			   xincco*(grid(6,n)-1)+Tinc*(grid(3,k)-1)+grid(2,j)
		      endif
		    else
		      id=Pinc*(grid(1,i)-1)+xincco*(idco2(grid(4,l)-1)-1)+xincco*idh2o(1,(grid(5,m)-1))+ &
		         xincco*(grid(6,n)-1)+Tinc*(grid(3,k)-1)+grid(2,j)
		    endif
	          else
		    id=Pinc*(grid(1,i)-1)+xincco*(idco2(grid(4,l)-1)-1)+xincco*idh2o((grid(4,l)-6),(grid(5,m)-1))+ &
		       xincco*(grid(6,n)-1)+Tinc*(grid(3,k)-1)+grid(2,j)
	          endif
	          read(kFile,rec=id) &
		    (kq(grid(1,i),grid(2,j),grid(3,k),grid(4,l),grid(5,m),grid(6,n))%values(iq),iq=1,Nq)
		endif
	      endif
	      !calculate k-distributions
	      final(1:Nq)=final(1:Nq)+weight*kq(grid(1,i),grid(2,j),grid(3,k), &
						grid(4,l),grid(5,m),grid(6,n))%values(1:Nq)
	    enddo !n
	  enddo !m
	enddo !l
      enddo !k
    enddo !j
  enddo !i

  return

end subroutine get_k_m

!-----------------------------------------------------------------------------------------

subroutine get_a(grid,Mix_Info,Trad,af)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Get a-values directly from XQtrTable and then carry on linear interpolation 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use commonData, only : GasMixInfo

  implicit none
  TYPE(GasMixInfo),INTENT(IN):: Mix_Info
  INTEGER,INTENT(IN) :: grid(6,num)
  real(8),intent(in) ::  Trad
  real(8),dimension(Nq),intent(out) :: af
  integer,parameter :: Pinc=274400
  integer :: i,j,k,l,m,n,id,iq

  af=0.d0
  wP=(Mix_Info%P-P(grid(1,1)))/(P(grid(1,2))-P(grid(1,1)))
  wTp=(Trad-Tp(grid(2,1)))/dTp
  wT=(Mix_Info%T-T(grid(3,1)))/dT
  if (Mix_Info%xCO2 < xCO2(6)) then
    wxco2=(Mix_Info%xCO2-xCO2(grid(4,1)))/dxco2_1
  else
    wxco2=(Mix_Info%xCO2-xCO2(grid(4,1)))/dxco2_2
  endif
  if (Mix_Info%xH2O < xH2O(6)) then
    wxh2o=(Mix_Info%xH2O-xH2O(grid(5,1)))/dxh2o_1
  else
    wxh2o=(Mix_Info%xH2O-xH2O(grid(5,1)))/dxh2o_2
  endif
  if (Mix_Info%xCO < xCO(2)) then
    wxco=(Mix_Info%xCO-xCO(grid(6,1)))/dxco_1
  elseif (Mix_Info%xCO < xCO(3)) then
    wxco=(Mix_Info%xCO-xCO(grid(6,1)))/dxco_2
  elseif (Mix_Info%xCO < xCO(4)) then
    wxco=(Mix_Info%xCO-xCO(grid(6,1)))/dxco_3
  else
    wxco=(Mix_Info%xCO-xCO(grid(6,1)))/dxco_4
  endif

  do i=1,num
    if (i==1) then; weP=1.d0-wP;
    else; weP=wP
    endif
    do j=1,num
      if (j==1) then; weTp=1.d0-wTp;
      else; weTp=wTp
      endif
      do k=1,num
        if (k==1) then; weT=1.d0-wT;
        else; weT=wT
        endif
	do l=1,num
          if (l==1) then; wexco2=1.d0-wxco2;
          else; wexco2=wxco2
          endif
	  do m=1,num
            if (m==1) then; wexh2o=1.d0-wxh2o;
            else; wexh2o=wxh2o
            endif
	    do n=1,num
              if (n==1) then; wexco=1.d0-wxco;
              else; wexco=wxco
              endif
	      !calculate weight
	      weight=weP*weTp*weT*wexco2*wexh2o*wexco
	      !load a-values into memory
	      if (associated(aq(grid(1,i),grid(2,j),grid(3,k), &
				grid(4,l),grid(5,m),grid(6,n))%values)) then
	        continue
	      else
	        allocate(aq(grid(1,i),grid(2,j),grid(3,k),grid(4,l),grid(5,m),grid(6,n))%values(Nq))
		id=Pinc*(grid(1,i)-1)+xincco2*(grid(4,l)-1)+xinch2o*(grid(5,m)-1)+ &
		   xincco*(grid(6,n)-1)+Tinc*(grid(3,k)-1)+grid(2,j)
	        read(aFile,rec=id) &
		  (aq(grid(1,i),grid(2,j),grid(3,k),grid(4,l),grid(5,m),grid(6,n))%values(iq),iq=1,Nq)
	      endif
	      !calculate a-values
	      af(1:Nq)=af(1:Nq)+weight*aq(grid(1,i),grid(2,j),grid(3,k), &
						grid(4,l),grid(5,m),grid(6,n))%values(1:Nq)
	    enddo !n
	  enddo !m
	enddo !l
      enddo !k
    enddo !j
  enddo !i

  return

end subroutine get_a

!-----------------------------------------------------------------------------------------

subroutine get_a_m(grid,Mix_Info,Trad,af)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!Get a-values directly from FullTable and then carry on linear interpolation 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use commonData, only : GasMixInfo

  implicit none
  TYPE(GasMixInfo),INTENT(IN):: Mix_Info
  INTEGER,INTENT(IN) :: grid(6,num)
  real(8),intent(in) ::  Trad
  real(8),dimension(Nq),intent(out) :: af
  integer,parameter :: Pinc=566048
  integer :: i,j,k,l,m,n,id,iq

  af=0.d0
  wP=(Mix_Info%P-P(grid(1,1)))/(P(grid(1,2))-P(grid(1,1)))
  wTp=(Trad-Tp(grid(2,1)))/dTp
  wT=(Mix_Info%T-T(grid(3,1)))/dT
  if (Mix_Info%xCO2 < xCO2(6)) then
    wxco2=(Mix_Info%xCO2-xCO2(grid(4,1)))/dxco2_1
  elseif (Mix_Info%xCO2 < xCO2(7)) then
    wxco2=(Mix_Info%xCO2-xCO2(grid(4,1)))/dxco2_2
  elseif (Mix_Info%xCO2 <= xCO2(nxco2)) then
    wxco2=(Mix_Info%xCO2-xCO2(grid(4,1)))/dxco2_3
  endif
  if (Mix_Info%xH2O < xH2O(6)) then
    wxh2o=(Mix_Info%xH2O-xH2O(grid(5,1)))/dxh2o_1
  elseif (Mix_Info%xH2O < xH2O(10)) then
    wxh2o=(Mix_Info%xH2O-xH2O(grid(5,1)))/dxh2o_2
  elseif (Mix_Info%xH2O <= xH2O(nxh2o)) then
    wxh2o=(Mix_Info%xH2O-xH2O(grid(5,1)))/dxh2o_3
  endif
  if (Mix_Info%xCO < xCO(2)) then
    wxco=(Mix_Info%xCO-xCO(grid(6,1)))/dxco_1
  elseif (Mix_Info%xCO < xCO(3)) then
    wxco=(Mix_Info%xCO-xCO(grid(6,1)))/dxco_2
  elseif (Mix_Info%xCO < xCO(4)) then
    wxco=(Mix_Info%xCO-xCO(grid(6,1)))/dxco_3
  elseif (Mix_Info%xCO < xCO(5)) then
    wxco=(Mix_Info%xCO-xCO(grid(6,1)))/dxco_4
  else
    wxco=(Mix_Info%xCO-xCO(grid(6,1)))/dxco_5
  endif

  do i=1,num
    if (i==1) then; weP=1.d0-wP;
    else; weP=wP
    endif
    do j=1,num
      if (j==1) then; weTp=1.d0-wTp;
      else; weTp=wTp
      endif
      do k=1,num
        if (k==1) then; weT=1.d0-wT;
        else; weT=wT
        endif
	do l=1,num
          if (l==1) then; wexco2=1.d0-wxco2;
          else; wexco2=wxco2
          endif
	  do m=1,num
            if (m==1) then; wexh2o=1.d0-wxh2o;
            else; wexh2o=wxh2o
            endif
	    do n=1,num
              if (n==1) then; wexco=1.d0-wxco;
              else; wexco=wxco
              endif
	      !calculate weight
	      weight=weP*weTp*weT*wexco2*wexh2o*wexco
	      !load a-values into memory
	      if (associated(aq(grid(1,i),grid(2,j),grid(3,k), &
				grid(4,l),grid(5,m),grid(6,n))%values)) then
	        continue
	      else
	        allocate(aq(grid(1,i),grid(2,j),grid(3,k),grid(4,l),grid(5,m),grid(6,n))%values(Nq))
		if (grid(4,l)<8) then
		  id=Pinc*(grid(1,i)-1)+xincco*idco2(grid(4,l)-1)+xincco*idh2o(1,(grid(5,m)-1))+ &
		     xincco*(grid(6,n)-1)+Tinc*(grid(3,k)-1)+grid(2,j)
		else
		  id=Pinc*(grid(1,i)-1)+xincco*idco2(grid(4,l)-1)+xincco*idh2o((grid(4,l)-6),(grid(5,m)-1))+ &
		     xincco*(grid(6,n)-1)+Tinc*(grid(3,k)-1)+grid(2,j)
		endif
	        read(aFile,rec=id) &
		  (aq(grid(1,i),grid(2,j),grid(3,k),grid(4,l),grid(5,m),grid(6,n))%values(iq),iq=1,Nq)
	      endif
	      !calculate a-values
	      af(1:Nq)=af(1:Nq)+weight*aq(grid(1,i),grid(2,j),grid(3,k), &
						grid(4,l),grid(5,m),grid(6,n))%values(1:Nq)
	    enddo !n
	  enddo !m
	enddo !l
      enddo !k
    enddo !j
  enddo !i

  return

end subroutine get_a_m

!-----------------------------------------------------------------------------------------

SUBROUTINE k_interp(nxy ,xx,yy, ni, xi,yi )

 IMPLICIT NONE
 INTEGER,INTENT(in) :: nxy, ni
 REAL(8),DIMENSION(nxy),INTENT(in) :: xx, yy
 REAL(8),DIMENSION(ni),INTENT(in) :: xi
 REAL(8),DIMENSION(ni),INTENT(out) :: yi

 INTEGER :: i, iq, ibgn
 REAL(8) :: dx

 ibgn=1
 DO iq= 1, ni
    
    ! constant extrapolate
    IF (xi(iq) < xx(1)) THEN
      yi(iq) = yy(1); cycle
    !ELSE IF (xi(iq) >= xx(nxy)) THEN
    !  yi(iq:ni)=yy(nxy); exit
    !END IF
    ! Not constant extrapolate
    elseif (xi(iq)>=xx(nxy)) then
      dx=xx(nxy)-xx(nxy-1)
      yi(iq)=yy(nxy)+(yy(nxy)-yy(nxy-1))*(xi(iq)-xx(nxy))/dx
      cycle
    endif

    DO i=ibgn, nxy
       IF (xx(i) > xi(iq) ) THEN ! this means xx(i-1)<=xi(iq)<xx(i)
      dx = xx(i)-xx(i-1)
          if (abs(dx)<tiny(xx)) then
        yi(iq)=yy(i-1)
      else
            yi(iq)= yy(i-1) + ( yy(i)-yy(i-1) ) * ( xi(iq)-xx(i-1) ) /dx
            ibgn=i
            EXIT
          endif
       END IF

    END DO ! i
 END DO ! iq

 RETURN
end subroutine

!-----------------------------------------------------------------------------------------

SUBROUTINE a_interp(nxy ,xx,yy, ni, xi,yi )

 IMPLICIT NONE
 INTEGER,INTENT(in) :: nxy, ni
 REAL(8),DIMENSION(nxy),INTENT(in) :: xx, yy
 REAL(8),DIMENSION(ni),INTENT(in) :: xi
 REAL(8),DIMENSION(ni),INTENT(out) :: yi

 INTEGER :: i, iq, ibgn
 REAL(8) :: dx

 ibgn=1
 DO iq= 1, ni
    
    ! constant extrapolate
    IF (xi(iq) < xx(1)) THEN
      yi(iq) = yy(1); cycle
    !ELSE IF (xi(iq) >= xx(nxy)) THEN
    !  yi(iq:ni)=yy(nxy); exit
    !END IF
    ! Not constant extrapolate
    elseif (xi(iq)>=xx(nxy)) then
      dx=xx(nxy)-xx(nxy-1)
      yi(iq)=yy(nxy)+(yy(nxy)-yy(nxy-1))*(xi(iq)-xx(nxy))/dx
      if (yi(iq)<0.d0) yi(iq)=yi(iq-1)
      cycle
    endif

    DO i=ibgn, nxy
       IF (xx(i) > xi(iq) ) THEN ! this means xx(i-1)<=xi(iq)<xx(i)
      dx = xx(i)-xx(i-1)
          if (abs(dx)<tiny(xx)) then
        yi(iq)=yy(i-1)
      else
            yi(iq)= yy(i-1) + ( yy(i)-yy(i-1) ) * ( xi(iq)-xx(i-1) ) /dx
            ibgn=i
            EXIT
          endif
       END IF

    END DO ! i
 END DO ! iq

 RETURN
end subroutine

!-------------------------------fsk_table-------------------------------------------

subroutine fsk_table (gasMix_f, refMix_f,Nq, k, a, g, gFSK, NoPDB, m)
! low level subroutine to calculate k and a vs g using FSK look-up table
use commonData, only : GasMixInfo


type(GasMixInfo),intent(in) :: gasMix_f,refMix_f !local and reference gas mixture information
integer, intent(in) :: Nq, NoPDB !number of quadrature points and constant number of points for FSK look-up table
integer, intent(in) :: m !table number
real(8), dimension(Nq), intent(inout) :: k, a, g !k and a values, quadrature points
real(8), dimension(NoPDB), intent(in) :: gFSK !quadrature points corresponding to data in FSK look-up table
!local variables
real(8), dimension(NoPDB) :: kpre, apre

if (m==1) then
  call kgaInterp(gasMix_f,refMix_f,refMix_f%T,kpre,apre,gFSK)
elseif (m==2) then
  call kgaInterpM(gasMix_f,refMix_f,refMix_f%T,kpre,apre,gFSK)
endif

call k_interp(NoPDB, gFSK, kpre,Nq, g, k)
call a_interp(NoPDB, gFSK, apre,Nq, g, a)

end subroutine fsk_table


end module fskTable
