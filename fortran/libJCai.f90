!!$**********************************************************************************
!!$                   SPECTRAL RADIATION CALCULATION SOFTWARE ~ SRCS ~
!!$**********************************************************************************
!!$
!!$                     COPYRIGHT (C) 2013 MICHAEL F. MODEST
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
module libCommon
! Procedures included in this module may already be included in other
! spectral modules. They are listed here to make libJCai stand-alone.
! Therefore, it is NOT suggest to directly include module "libCommon"
! with a "use" statement unless not other spectral modules are present.
! It is wise to control access using "use ..., only:..." statement.

implicit none
! Table of Contents
public :: &
 BBFN,& ! Fractional blackbody emission power function 
 quadgen2 ! Quadrature generator
         

contains
!-----------------------------------------------------------------------
! black body emission function from absc_kg_gen.f90
  REAL(8) FUNCTION BBFN(Q) 
!!$     ********************************************************************
!!$     *  This subroutine calculates the fractional blackbody             *
!!$     *  emissive power f(n*lambda*T), where X=n*lambda*T in (micro-m*K) *
!!$     ********************************************************************
    REAL(8) :: PI,CC,C2,EPS,V,EX,M,EM,VM,BM,Q
    REAL(8),PARAMETER :: x1= 21.d0, x2= 75671.d0

    IF (Q <x1) THEN        ! if Q is too small, return 0.0
       BBFN=0.d0; RETURN
    ELSEIF (Q>x2) THEN     ! if Q is too large, return 1.0
       BBFN= 1.d0; RETURN
    ENDIF   

    PI=4.D0*DATAN(1.D0)
    CC=1.5D1/PI**4
    C2=1.4388D4
    EPS=1.D-16

    V=C2/Q  
    EX=DEXP(V)

    M=0
    BBFN=0.D0
    EM=1.D0 
5   M=M+1   
    VM=M*V  
    BM=(6.D0+VM*(6.D0+VM*(3.D0+VM)))/M**4
    EM=EM/EX
    BBFN=BBFN+BM*EM
    IF(VM**3*EM.GT.EPS) GOTO 5
    BBFN=CC*BBFN
    RETURN  
  END FUNCTION BBFN
!------------------------------------------------------------------------
! quadrature generator from quadlib.f90
SUBROUTINE quadgen2(Cheb2,g,w,NoP,alpha) ! quadrature open at both ends
  IMPLICIT NONE
  LOGICAL,INTENT(IN)    :: Cheb2
  INTEGER,INTENT(IN)    :: NoP  
  REAL(8),INTENT(IN),optional    :: alpha  ! transformation factor
  REAL(8),INTENT(OUT)   :: g(NoP),w(NoP)
!local variables
  REAL(8) :: gg(2*NoP),ww(2*NoP),dummy
  INTEGER(2) :: i

    CALL gaucheb2(gg,ww,2*NoP)
    g= -gg(NoP+1:2*NoP)
    w= ww(NoP+1:2*NoP)
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
END SUBROUTINE quadgen2
!------------------------------------------------------------------
! Gaussian quadrature calculation subroutine from quadlib.f90
SUBROUTINE gaucheb2(x,w,n)
IMPLICIT NONE
  INTEGER :: n
  DOUBLE PRECISION,PARAMETER :: PI= 3.1415926535897932384626
  REAL(8) :: x(n),w(n)
  DOUBLE PRECISION :: theta, sum
  INTEGER :: m,k
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx> By: Anquan Wang <xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
! DISCRIPTION:                                                                                       !
!   Calculate the zeros and weights of Tchebysheff polynomial of 2nd kind                            !
!   The zeros can be used as abscissas of Gaussian quadrature which are distributed between (-1,1)   !
!   The integral can be evaluated in the form:                                                       !
!                         Integral(f(x)*dx)= Sum(Wi*f(xi))                                           !
!----------------------------------------------------------------------------------------------------!
! PARAMETERS:                                                                                        !
!   n    -> the degree of Tchebyscheff polynomial of 2nd kind                                        !
!   x(n) -> the array of abscissas                                                                   !
!   w(n) -> the array of weights                                                                     !
!----------------------------------------------------------------------------------------------------!
! MEMO:                                                                                              !
!   P. J. Davis and P. Rabinowitz, Methods of Numerical Integration, 2nd edt, pp85, eqn(2.5.5.6,.8)  !
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx> 2003-06-02 <xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
  DO k= 1,(n+1)/2
    theta= REAL(k)*PI/(n+1.d0)
    x(k)= DCOS(theta)
    sum= 0.d0
    DO m= 1,(n+1)/2
      sum= sum+DSIN((2*m-1)*theta)/(2.d0*m-1.d0)
    ENDDO   
    w(k)= 4.d0*DSIN(theta)*sum/(n+1.d0)
    x(n+1-k)= -x(k)
    w(n+1-k)= w(k)
  ENDDO
END SUBROUTINE gaucheb2

end module libCommon

!*****************************************************************************
! new spectral utility functions
module libJCai
! Utility functions for spectral calculation
implicit none
integer,parameter,private ::dp=kind(1.d0)
! Table of Contents
Public :: &
 kPowerLaw,& ! generate a list of k values between kmin and kmax from power-law
 linearInterpMono,& ! linear interpolation 
 multiplicMix,& ! mix k-dist with multiplication rule
 MRmix,& ! mix two k-dist with Modest&Riazzi rule
 nbk2fsk,&   ! narrowband k-dist to full-spectrum k-dist
 superposiMix,& ! mix k-dist with superposition rule
 fsk_corr,&  ! fsck from correlation table of CO2 and H2O
 afun  ! a new "a" function calculation routine
contains
function kDist (absc, eta, w, k) result (g)
! low level general k-distribution calculation function 
! It is assumed that integratioin of w over eta is unity.
! For discrete data integration of w over eta means
!    sum_i w_i*(eta(i+1)-eta(i))
! Variable eta is assumed to consecutive and increasing.
! If the integration of w over eta is less than unity, then g
real(dp), dimension(:), intent(in) :: absc
real(dp), dimension(size(absc,1)), intent(in) :: eta
real(dp), dimension(size(absc,1)-1), intent(in) :: w
real(dp), intent(in) :: k ! scalar value k
real(dp) :: g
! local variables
real(dp) :: frac 
  ! w weighted fraction of spectral interval
!real(dp), dimension(size(absc,1)-1) :: dEta ! spectral interval width
logical, dimension(size(absc, 1)) :: isGT
integer :: n, i
n = size(absc,1)
do i = 1, n
isGT(i) = k .GT. absc(i)
enddo
!dEta = eta(2:n)-eta(1:(n-1))
! count the fraction of spectral interval with absc>k
g = 0.d0
do i=1, n-1
  if (isGT(i)) then ! k>absc(i)
    if (isGT(i+1)) then ! k>absc(i) & k>absc(i+1): interval does not count
      frac = 0.d0
    else  ! k>absc(i) & k<=absc(i+1): interval partially counts
      frac = (absc(i+1)-k)/(absc(i+1)-absc(i))*(eta(i+1)-eta(i))*w(i)
    endif
  else  ! k<=absc(i)
    if (isGT(i+1)) then ! k<=absc(i) & k>absc(i+1): interval partially counts
      frac = (absc(i)-k)/(absc(i)-absc(i+1))*(eta(i+1)-eta(i))*w(i)
    else ! k<=absc(i) & k<=absc(i+1): interval fully counts
      frac = (eta(i+1)-eta(i))*w(i)
    endif
  endif
  g = g + frac
enddo
g = 1.d0 - g 
! Notes on the result
! (1) if k>maxval(abac), g = 1, i.e. g never exceeds 1 for non-negative w and 
!     increasing eta
! (2) if k<minval(absc), g = 1-sum(w*dEta) 
! (3) If sum(w*dEta)<1 (e.g. with spectral leakage), g>gmin for k>0, where
!     gmin = 1-sum(w*dEta). This is equivalent to assuming absc=0 over 
!     out-of-range spectrum
! (4) If sum(w*dEta)>1 which is unphysical, g<0 for small k, so the output
!     is also unphysical to warn user. 
! (5) There is no intension for this function to do any special treatments
!     such as bounding g between 0 and 1. Output is for error probing.
end function kDist
!--------------------------------------------------------------
function kPowerLaw (kmin,kmax, n, pwr) result(k)
! function generate a list of k values between kmin and kmax
! according to power law with power "pwr"
real(dp),  intent(in):: kmin, kmax, pwr
integer, intent(in) :: n
real(dp), dimension(n) :: k
real(dp) :: pwrk_min, pwrk_max, pwrk_step
integer :: i
pwrk_min = kmin**pwr
pwrk_max = kmax**pwr
pwrk_step = (pwrk_max-pwrk_min)/real(n-1, dp)
k = (/(pwrk_min+real(i-1,dp)*pwrk_step, i=1, n)/)
k = k**(1.d0/pwr)
end function kPowerLaw

!------------------------------------------------------------
function LBL2FSK_k (absc, eta, T, k) result(g)
! Calculate full-spectrum k-distribution from line-by-line absorption 
! coefficient for a fixed list of k values
use libCommon, only : BBFN
real(dp), intent(in), dimension(:) :: absc, k
! absorption coefficient, k values for g
real(dp), intent(in), dimension(size(absc,1)) :: eta ! wavenumber
real(dp), intent(in) :: T ! temperature
real(dp), dimension(size(k,1)) :: g
real(dp), dimension(size(eta,1)) :: Q, BB 
  ! dimensionless spectral varible, black body emission power 
  ! See BBFN and Appendix C of Modest 2003 for more details 
real(dp), dimension(size(absc,1)-1) :: Ib
  ! Planck function calculated by Delta Q/Delta eta
  ! Will give the reason in doc
integer :: i, nk, n
nk = size(k,1)
n = size(eta,1)
Q = T*1.d4/eta
do i = 1, n ! loop over all wavenumbers for black body emission power
  BB(i) = BBFN(Q(i))
enddo 
! Original function calls LBL2FSK_sk for each k
! It turns out that calculation of BB takes the most of time
! So expand LBL2FSK_sk inline to avoid duplicated calculation of BB
Ib = (BB(1:(n-1))-BB(2:n))/(eta(2:n)-eta(1:(n-1)))
do i = 1, nk
  g(i) = kDist (absc, eta, Ib, k(i))
!  print*, i
enddo
end function LBL2FSK_k
!-------------------------------------------------------------
function LBL2FSK_n (absc, eta, T, n) result(gk)
! Calculate full-spectrum k-distribution from line-by-line absorption 
! coefficient for over n k values, the list of k is distributed from
! power law between minimum and maximum of absorption coefficients 
real(dp), intent(in), dimension(:) :: absc
! absorption coefficient
real(dp), intent(in), dimension(size(absc,1)) :: eta ! wavenumber
real(dp), intent(in) :: T ! temperature
integer, intent(in) :: n ! total number of k's
real(dp), dimension(n,2) :: gk
! output is two column vectors, first column for g, second for k
! in ascending order
!real(dp), dimension(size(eta,1)) :: Q, BB 
  ! dimensionless spectral varible, black body emission power 
  ! See BBFN and Appendix C of Modest 2003 for more details 
!real(dp), dimension(size(absc,1)-1) :: Ib
  ! Planck function calculated by Delta Q/Delta eta
  ! Will give the reason in doc
real(dp) :: kmin, kmax
real(dp), parameter :: pwr = 0.1d0, kminBound = 1.d-9
kmin = minval(absc,1); kmax = maxval(absc,1);
kmin = max(kmin, kminBound) ! bound the mininum k
gk(:,2) = kPowerLaw (kmin,kmax, n, pwr)
gk(:,1) = LBL2FSK_k(absc, eta, T, gk(:,2))
end function LBL2FSK_n
!------------------------------------------------------------
function LBL2FSK_sk (absc,eta, T, k) result(g)
! Calculate full-spectrum k-distribution from line-by-line absorption
! coefficient at a single k value
! This is low level calculation routine for 
use libCommon, only : BBFN
real(dp), intent(in), dimension(:) :: absc
real(dp), intent(in), dimension(size(absc,1)) :: eta
real(dp), intent(in) :: T, k ! temperature, scalar k
real(dp) :: g
! local variables
real(dp), dimension(size(absc,1)) :: Q, BB 
  ! dimensionless spectral varible, black body emission power 
  ! See BBFN and Appendix C of Modest 2003 for more details 
real(dp), dimension(size(absc,1)-1) :: Ib
  ! Planck function calculated by Delta Q/Delta eta
  ! Will give the reason in doc
integer :: i, n
n = size(eta,1)
Q = T*1.d4/eta
do i = 1, n ! loop over all wavenumbers for black body emission power
  BB(i) = BBFN(Q(i))
enddo
Ib = (BB(1:(n-1))-BB(2:n))/(eta(2:n)-eta(1:(n-1)))
g = kDist (absc, eta, Ib, k)
end function LBL2FSK_sk
!------------------------------------------------------------
function LBL2NBK_fs (absc, eta, interval, n) result(gk)
! Calculate narrowband k-dist from line-by-line database for a list of
! spectral invervals. Inside of this program, the validity of narrowband
! assumption (i.e. small variation of the Planck function) is not checked.
implicit none
real(dp), intent(in), dimension(:) :: absc ! absorption coefficient
real(dp), intent(in), dimension(size(absc,1)) :: eta ! wavenumber
real(dp), intent(in), dimension(:, :) :: interval 
  ! list of intervals, colume 1 is lower bound, colume 2 is upper bound
integer, intent(in) :: n
real(dp), dimension(size(interval,1), n, 2) :: gk
  ! output for narrowband k-dist, the first index is for narrowband, second
  ! for the total number of k values, third index is 1 for g, 2 for k.
! local variables
integer :: i, ik, jn, nnb, inb, nEta, lb, ub
real(dp) :: etaMin, kmax, kmin
real(dp), dimension(2) :: abscInt ! absc at interval points
!real(dp), dimension(size(interval)*2) :: interval1, absc1
  ! temporary array for interval bounds 
real(dp), dimension(:), allocatable :: nbEta, nbAbsc, nbW
logical, dimension(size(interval,1),2) :: mask
real(dp), parameter :: eps = 1.d1*tiny(1.d0) ! give a small number for nonexact real
real(dp), parameter :: pwr = 0.1d0, kminBound = 1.d-9
nnb = size(interval,1) 
gk = 0.d0
mask =  interval.lt.(eta(1)-eps)
if (any(any(mask,2),1)) then
  write(*,*) "Lowest wavenumber is not capable to cover all intervals"
  stop
end if
mask = interval.gt.(eta(size(eta,1)) +eps)
if (any(any(mask,2),1)) then
  write(*,*) "Greatest wavenumber is not capable to cover all intervals"
  stop
end if
!inverval1 = reshape(interval(1:nnb,1:2), (/nnb*2/))
nEta = size(eta,1)
!call linearInterpMono(nEta, eta, absc, nnb*2, interval1, absc1)
!abscInt = reshape(absc1, (/nnb, 2/))
do inb = 1, nnb 
  ! identify spectral segment
  lb=1; ub=neta
  do i = 1, neta
    if (eta(i).ge.interval(inb,1)) then
      lb = i;exit
    endif
  enddo
  do i = lb, neta
    if (eta(i).ge.interval(inb,2)) then
      ub = i-1; exit
    endif
  enddo
  if (ub.le.lb) then 
    write(*,*) "Interval ", inb, " lower bound ", interval(inb,1), " and upper bound ", &
      interval(inb,2), " are with in the same spectral interval" 
    stop
  endif
  call linearInterpMono(neta, eta, absc, 2, interval(inb,1:2), abscInt)
  allocate(nbEta(ub-lb+3), nbAbsc(ub-lb+3), nbw(ub-lb+2))
  nbEta = (/interval(inb,1),eta(lb:ub), interval(inb,2)/)
  nbAbsc = (/abscInt(1), absc(lb:ub), abscInt(2)/)
  nbW = 1.d0/(interval(inb, 2)-interval(inb,1))
  kmax=maxval(nbAbsc); kmin=minval(nbAbsc)
  kmin = max(kmin, kminBound) ! bound the mininum k
  kmax = max(kmax, kmin + kminBound) ! bound the max to be at least 2*kminBound, and larger than kmin
  gk(inb,:,2) = kPowerLaw (kmin,kmax, n, pwr)
  do ik = 1, n
    gk(inb, ik, 1) = kDist(nbAbsc, nbEta, nbW, gk(inb, ik,2))
  end do
  deallocate(nbEta, nbAbsc, nbW) 
end do
end function LBL2NBK_fs
!-----------------------------------------------------------------------
subroutine linearInterpMono(nxy, xx, yy, ni, xi, yi) 
! simple linear interpolation using closest points
! constant value for extrapolation
! assume xx and xi are monotonically increasing
! (1) It is assumed that xx is in ascendant order, which is typical for
!     linear interpolation
! (2) xi is also assumed to be in ascendant order, which improves search
!     speed
! (3) Constant extrapolation using two boundary values.
integer, intent(in) :: nxy, ni 
integer, parameter :: dp = kind(1.d0)
real(dp), dimension(nxy), intent(in) :: xx, yy
real(dp), dimension(ni), intent(in) :: xi
real(dp), dimension(ni), intent(out) :: yi
integer :: i, n 
n = 1
loopi:do i = 1, ni 
  do while (xi(i).ge.xx(n))
    ! xi(i) is not in interval xx(n-1) to xx(n)
    n = n + 1 ! then move to next interval
    if (n .gt. nxy) then ! xi(i) is larger than the largest xx
      ! out of bound and use constant value
      yi(i:ni) = yy(nxy) 
      exit loopi ! no need to perform further calculation
    end if  
  end do  
  ! now xi(i) is between xx(n-1) and xx(n), except xi(i) is smaller than xx(1)
  if (n.eq.1) then 
    ! xi(i) is smaller than the smallest of xx
    yi(i) = yy(1) 
  else 
    yi(i) = yy(n)+(yy(n-1)-yy(n))*(xi(i)-xx(n))/(xx(n-1)-xx(n))
  end if  
end do loopi
end subroutine linearInterpMono

!------------------------------------------------------------
function multiplicMix(g) result(gmix)
! Mix k-dist with the multiplication rule
real(dp), dimension(:,:), intent(in) :: g
real(dp), dimension(size(g,1)) :: gmix 
gmix = product(g, 2)
end function multiplicMix

!------------------------------------------------------------
function MRmix(k1, g1, k2, g2, k, quadNop, quadScheme) result(g)
use libCommon, only: quadgen2
real(dp), dimension(:), intent(in) :: k1, k2, k
real(dp), dimension(size(k1,1)), intent(in) :: g1
real(dp), dimension(size(k2,1)), intent(in) :: g2
integer, intent(in), optional :: quadScheme, quadNop
  ! user-chosen scheme and number-of-points for integration
real(dp), dimension(size(k, 1)) :: g
real(dp), dimension(:), allocatable :: gq1,gq2, wq1 
  ! quad points and weights to perform numerical integration
real(dp), dimension(:), allocatable :: kq1, kq2 
  ! k1 and k2 interpolate at quad points
integer :: i, nk1, nk2, nop, scheme
real(dp) :: kmax, kmin
integer, parameter :: defaultQuadNop=64, defaultQuadScheme=0
! (1) note that integration is over g1, so if g1 is defined over quadrature points
! then interpolation has no errors
! (2) If multiple species are mixed and they are databased on quadrature points
! then it is possible to mix them one by one, and input previously mixed k-g as k2 g2
! (3) There are different ways to perform integration numerically, user may choose
! one by specifying "scheme" and "order"
! (4) k1 and k2 can be of different sizes
! (5) if other integration rules are needed, just add a new scheme number and a different
! set of weights
! (6) By default 64 points quadrature open at both ends will be used for integration
if (present(quadNop)) then
  nop = quadNop
else 
  nop = defaultQuadNop
endif
if (present(quadScheme)) then
  scheme = quadScheme
else 
  scheme = defaultQuadScheme
endif
allocate(gq1(nop),gq2(nop),wq1(nop),kq1(nop),kq2(nop))
select case (scheme)
  case (0)
    call quadgen2(.false., gq1, wq1, nop)
  case default
    print*, 'Not implemented yet'; stop
end select
g = 1.d0 ! by default let g always be one 
nk1 = size(k1,1); nk2 = size(k2,1)
! interpolate g1 to quad points
call linearInterpMono  (nk1, g1, k1, nop, gq1, kq1)
kmin = k1(1) + k2(1) ! usually k1 and k2 are monotonic
kmax = k1(nk1) + k2(nk2)

do i = 1, size(k,1)
  if (k(i) .lt. kmin) then
    g(i) = 0.d0; cycle    ! k is too small, next loop
  endif
  if (k(i) .gt. kmax) exit ! k is too large, no need to loop other k's
  kq2 = k(i) - kq1
  ! kq2 is in descent order
  call linearInterpMono  (nk2, k2, g2, nop, kq2(nop:1:-1), gq2(nop:1:-1))
  gq2 =max(gq2, 0.d0) ! lower bound zero
  g(i) = sum(gq2*wq1)
enddo
deallocate(gq1, gq2, wq1, kq1, kq2)
end function MRmix

!-------------------------------------------------------------
function nbk2fsk(T, gnb, intnb) result (g)
! function converts narrowband k-distribution into full-spectrum
! k-distribution at given temperature T.
! Input:
! T: double precision scalar, temperature
! gnb: double precision NxL matrix, narrowband k-distributions.
!    : (1) The first index is for narrowband id, the second index
!    : is for k-distribution data. So gnb has N narrowbands and
!    : each narrowband has L datapoints for k-distribution.
!    : (2) It is assumed that all narrowband g values are defined
!    : on the same k values, this is why k values are not 
!    : required in calculation.
!    : (3) Note that g values are not required to be monotonic. 
!    : which means the k values that g is defined on can be in 
!    : arbitrary order. This is because the calculation at each
!    : fixed k is indepedent of other k values. But all narrowband
!    : g values must be of the same order with k. In program this 
!    : means
!    :      g(jN, jL) is the g value at k(jL) for any jN, jL
!    : (4) Although in real applications k as well as g are all in
!    : non-decreasing or non-increasing order, such structure is
!    : not used
!    : (5) The narrowband can be listed in any order, as long as
!    : being listed in the same order as narrowband intervals
! intnb: double precision Nx2 matrix, narrowband interval bounds.
!    : (1) The list of intervals must be given in the same order
!    : as gnb, in program this means
!    :      g(jN,:) is narrowband k-g over spectral interval 
!    :      [intnb(jN,1), intnb(jN,2)] for any jN
!    : (2) The first column must be interval lower bounds, while
!    : the second column must be interval upper bounds. However,
!    : this order is not checked in program. Reversing this order
!    : will cause negative g values.
!    : (3) The major reason to develop this function is handling
!    : spectral leakage due to limited spectral ranges. 
!    : (4) Since spectral variable wavenumber shows in denominator
!    : in dimensionless spectral variable for function BBFN, 
!    : lowest wavenumber should not be zero, but a small number.
use libCommon, only : BBFN
implicit none
integer, parameter :: dp = kind(1.d0)
real(dp), intent(in) :: T ! temperature
real(dp), dimension(:,:), intent(in) :: gnb ! narrowband g
real(dp), dimension(size(gnb,1), 2), intent(in) :: intnb 
    ! narrowband interval bounds, first column must be lower
    ! bounds, second column must be upper bounds
real(dp), dimension(size(gnb,2)) :: g 
    ! full spectrum g as output

! local variables
real(dp), dimension(size(gnb,1), 2) :: bb 
    ! dimensionless spectral variable for Planck function
real(dp), dimension(size(gnb,1)) :: w 
    ! spectral weight from Planck function
integer :: jN, jL, N, L
real(dp), parameter :: maxLeak=1.d-2 
    ! maximum allowed leakage is 1%
g = 0.d0; bb = (T*1.d4)/intnb ! initialize
N = size(gnb,1); L = size(gnb,2) ! prepare loop bounds
do jN = 1, N
  w(jN) = BBFN(bb(jN,1))-BBFN(bb(jN,2)) ! find spectral weight
enddo
if (abs(sum(w)-1.d0)>maxLeak) & 
  write(*,*) "Sum of spectral weight deviates more than 1% from 1."
! give a warning for possible spectral leakage
do jL = 1, L
  g(jL) = sum(w*(1.d0-gnb(:,jL))) 
  ! count larger k's, so over leakage spectrum k is assumed to be 0
end do
g = 1.d0-g
end function NBK2FSK

!--------------------------------------------------------------
function superposiMix(g) result(gmix)
! Mix k-dist with the superposition rule
real(dp), dimension(:,:), intent(in) :: g
real(dp), dimension(size(g,1)) :: gmix 
integer :: m
m = size(g,2)
gmix = sum(g, 2)
gmix = gmix +1.d0 - real(m, dp)
end function superposiMix


SUBROUTINE Fskdh2oms(Tg,Tb,k,x,gcal)
   ! modest & singh correlation for H2O
   IMPLICIT NONE
   INTEGER  :: l,m,n
   DOUBLE PRECISION :: bigp,xi,xis,gcal,Tg,Tb,k,x
   DOUBLE PRECISION :: almn(0:3,0:3,0:3),blmn(0:1,0:2,0:2)
   DATA almn /   0.11727390382594D+01,  -0.79233701090528D+00, &
         0.45678148180508D+00,  -0.67355817198147D-01, &
         0.88261090723411D+00,   0.33022055332160D+01, &
        -0.20480774079693D+01,   0.39413787303095D+00, &
        -0.28464957949024D+00,  -0.16339745840812D+01, &
         0.10992463891114D+01,  -0.20874526756693D+00, &
         0.52641314250448D-01,   0.26822395746801D+00, &
        -0.18705727454161D+00,   0.35072746249398D-01, &
         0.27995556315637D+00,   0.52705545762030D+00, &
        -0.43155723889914D+00,   0.13503443577078D+00, &
         0.13110405345336D+01,  -0.13535195345260D-01, &
         0.23946512473986D-01,  -0.37574344913040D-01, &
        -0.61658230571746D+00,  -0.81576227057821D-01, &
         0.15935609591997D+00,  -0.22954767376049D-01, &
         0.10234050877759D+00,   0.17012916696670D-01, &
        -0.38405881155705D-01,   0.65175870635333D-02, &
         0.86277669711840D-01,   0.34948997482539D+00, &
        -0.20189360135442D+00,   0.55391027176244D-01, &
         0.31365825397722D+00,  -0.10602534655868D+01, &
         0.85700017183181D+00,  -0.17833034900819D+00, &
        -0.13677988714456D+00,   0.57680694002893D+00, &
        -0.48058859578561D+00,   0.91345067362557D-01, &
         0.20072348102088D-01,  -0.97271482624167D-01, &
         0.81344362337134D-01,  -0.14501502880846D-01, &
         0.55203164526025D-01,  -0.90309091713057D-01, &
         0.22803582483225D+00,  -0.50412939547401D-01, &
        -0.26467885332271D-01,  -0.63565757156191D-01, &
        -0.13107199511879D+00,   0.37171590780014D-01, &
         0.13331890848674D-01,   0.67414123184757D-01, &
         0.25816812913826D-01,  -0.13250915677532D-01, &
        -0.25024127061198D-02,  -0.13749665053479D-01, &
        -0.14234808917787D-02,   0.19916213085042D-02/
   
   DATA blmn /  -0.89871032164934D+00,   0.53900322055117D+00, &
         0.10116065458653D+01,  -0.58957492617859D+00, &
        -0.48279630261745D+00,   0.27472330190048D+00, &
         0.10985053030604D+01,  -0.58298454159600D-01, &
        -0.10454377949888D+01,   0.11542729243720D+00, &
         0.47174494979936D+00,  -0.56228311600148D-01, &
        -0.22529697071712D+00,  -0.60603363774397D-02, &
         0.18432229208060D+00,  -0.16727319474944D-01, &
        -0.76913213683316D-01,   0.10350069781625D-01/
 

   xi=LOG10(k)  ! log 10
   xis=0d0  ! if x=0, xis=0 because of x**(m+1)
   DO l=0,2
     DO m=0,2
       DO n=0,1
         xis=xis+blmn(n,m,l)*(Tg/1000d0)**l*xi**n*x**(m+1)
       ENDDO !n
     ENDDO !m
   ENDDO !l
   bigp=0d0
   DO l=0,3
     DO m=0,3
       DO n=0,3
         bigp=bigp+almn(n,m,l)*(Tg/1000d0)**n*(Tb/1000d0)**m*(xi+xis)**l
       ENDDO !n
     ENDDO !m
   ENDDO !l
   gcal=0.5d0*tanh(bigp)+0.5d0
   RETURN
   END SUBROUTINE Fskdh2oms


   SUBROUTINE Fskdco2mm(Tg,Tb,absco,gcal)
   ! modest & mehta correlation for CO2
   IMPLICIT NONE
   INTEGER  :: l,m,n
   DOUBLE PRECISION :: bigp,ksai,gcal,Tg,Tb,absco
   DOUBLE PRECISION :: dlmn(0:3,0:3,0:3)
  DATA dlmn /  1.85071   ,  0.33373   ,  0.62660   , -0.12890   , & ! n=0~3, m=0, l=0
           -0.20643   , -2.57690   ,  0.30090   , -0.14090   , & ! n=0~3, m=1, l=0
            0.27664   ,  1.81420   , -0.24728   ,  0.10052   , & ! n=0~3, m=2, l=0
           -0.37435d-1, -0.37762   ,  0.53014d-1, -0.20836d-1, & ! n=0~3, m=3, l=0
            0.67523   ,  1.25760   ,  0.67523d-1, -0.39669d-1, & ! n=0~3, m=0, l=1
           -0.70897   , -3.07080   ,  1.71150   , -0.57694   , & ! n=0~3, m=1, l=1
            0.48493   ,  2.04603   , -1.20220   ,  0.40510   , & ! n=0~3, m=2, l=1
           -0.98138d-1, -0.41928   ,  0.24956   , -0.84109d-1, & ! n=0~3, m=3, l=1
            0.20690   ,  0.28500   , -0.48324d-1,  0.15174d-1, & ! n=0~3, m=0, l=2
           -0.39473   , -0.42333   ,  0.54095   , -0.23469   , & ! n=0~3, m=1, l=2
            0.29020   ,  0.21882   , -0.34748   ,  0.15857   , & ! n=0~3, m=2, l=2
           -0.61998d-1, -0.38629d-1,  0.68728d-1, -0.32314d-1, & ! n=0~3, m=3, l=2
            0.38488d-1,  0.18292d-1, -0.18958d-1,  0.61307d-2, & ! n=0~3, m=0, l=3
           -0.41013d-1, -0.59115d-2,  0.47118d-1, -0.23229d-1, & ! n=0~3, m=1, l=3
            0.37740d-1, -0.30114d-1, -0.11303d-1,  0.12140d-1, & ! n=0~3, m=2, l=3
           -0.87906d-2,  0.98357d-2, -0.53884d-4, -0.20604d-2/
   ksai=LOG10(absco)
   bigp=0d0
   DO l=0,3
     DO m=0,3
       DO n=0,3
         bigp=bigp+dlmn(n, m, l)*(Tg/1000d0)**n*(Tb/1000d0)**m*ksai**l
       ENDDO !n
     ENDDO !m
   ENDDO !l
   gcal=0.5d0*tanh(bigp)+0.5d0
   RETURN
   END SUBROUTINE Fskdco2mm


subroutine fsk_corr (Tg, xCO2, xH2O, Tref, xCO2ref, xH2Oref,& 
            & mixModel, Nq, g,w, k, a, errflag) 
! Top user level routine for 
! full spectrum correlated k-distribution including only H2O and CO2 using
! correlation tables
! WARNING: pressure must be 1 bar
! Notes: need only two temperatures "current" and "reference"
! How to map planck function temperature? Tp = Tg or Tref
use libCommon, only: quadgen2
integer, intent(in) :: mixModel, Nq
real(dp), intent(in) :: Tg, Tref, xCO2, xH2O, xCO2ref, xH2Oref
real(dp), dimension(nq), intent(in) :: g,w
real(dp), dimension(Nq), intent(out) :: k, a
integer, intent(out) :: errflag
! local variables
integer :: nq2, ik, mix! internal quadrature numbers for individual species before mixing
real(dp), dimension(:), allocatable :: gcr, grc, grr,gcc,k2, af, afs
real(dp), dimension(size(g,1)) :: kr, gr!, gq,wq
! TODO kcr, krc, krr needed?
! three k-dist needed for intermediate calculations
! kcr, gcr: k,g with absc eval at current thermo state, weighted by Planck fun at ref temperature
! krc, grc: k,g with absc eval at reference thermo state, weighted by Planck fun at current temperature
! krr, grr: k,g with absc eval at ref thermo state, weighted by Planck fun at ref temperature
real(dp), parameter :: kmin=1.d-9, kmax=1d4, pwr=0.1d0, smallNumber=1d-14, smoothPara=1d-6
integer, parameter :: expanRatio = 4 ! this expansion ratio controls accuracy if interpolation is needed.
! TODO data validation checking
mix = mixModel
! call quadgen2(.true., g, w, Nq)
! TODO need to remove quad 
a=0.d0 ! wired-in a functions
nq2 = min (128, Nq*expanRatio)
! internal calculation involves much longer k-dist to retain accuracy after interpolation.
allocate(gcr(nq2),grr(nq2),grc(nq2), k2(nq2), af(nq2), afs(nq2), gcc(nq2))
k2=0.d0; grc=0.d0; grr=0.d0; gcr=0.d0;af=0.d0; afs=0.d0
k2=kPowerLaw(kmin, kmax, nq2, pwr)
! three k-dist are involved in total
! the first one k-dist evaluated at current thermo state, with planck function at Tref
! Result gives k^* from correlated assumption
!call fsk_corr_mix(xCO2, xH2O, Tg, Tref, k2,gcr, mix, errflag)

! second one k-dist evaluated at ref thermo state, with planck function at current T
call fsk_corr_mix(xCO2ref, xH2Oref, Tref, Tg, k2,grc, mix, errflag)
! third one k-dist evaluated at ref thermo state, with planck function at Tref
call fsk_corr_mix(xCO2ref, xH2Oref, Tref, Tref, k2,grr, mix, errflag)

! calculate "a" function from ratio of above two k-dist
!! old a function
!af(1) = (grc(2)-grc(1))/(grr(2)-grr(1)+smallNumber)
!af(nq2) = (grc(nq2)-grc(nq2-1))/(grr(nq2)-grr(nq2-1)+smallNumber)
!do ik = 2, nq2-1
!    af(ik) = (grc(ik+1)-grc(ik-1))/(grr(ik+1)-grr(ik-1)+smallNumber)
!enddo
!call smooth_spline(nq2, nq2, k2, af, smoothPara,afs) 
!call linearInterpMono(nq2, grr, afs, Nq, g, a)
!! new a fucntion
!call quadgen2(.true., gq, wq, nq)
call afun(grr,k2, grc,k2, g, w,a)
! call interpolation to the required values
! old k^*
!call linearInterpMono(nq2, gcr, k2, Nq, g, k)
! new k^*
call linearInterpMono(nq2, grr, k2, nq, g, kr)
call linearInterpMono(nq2, k2, grc, nq, kr, gr)
call fsk_corr_mix(xCO2, xH2O, Tg, Tg, k2, gcc, mix, errflag)
call linearInterpMono(nq2, gcc, k2, nq, gr, k)

deallocate(k2, gcr,grc,grr,af, afs)
end subroutine fsk_corr

subroutine fsk_corr_mix(xCO2, xH2O, Tg, Tp, k, g, m, errflag)
! low level subroutine to calculate k vs g for a mixture with mole fraction
! xCO2, xH2O and Tg for gas temperature, Tp gives Planck function temperature
real(dp), intent(in) :: xCO2, xH2O, Tg, Tp
integer, intent(in) :: m  ! mixing model 
integer, intent(out) :: errflag ! error flag
real(dp), dimension(:), intent(in) :: k
real(dp), dimension(size(k,1)), intent(out) :: g
! local variables
integer :: nq, ik, mix
real(dp), dimension(size(k,1)) :: gCO2, gH2O
real(dp) :: kp ! pressure based k value, not planck mean
nq = size(k,1)
! TODO double array sizes for MRMix?
mix=m;gCO2=1.d0;gH2O=1.d0
if (xCO2.gt.0.d0) then
    do ik=1, nq
        kp=k(ik)/xCO2
        call fskdco2mm(Tg, Tp, kp, gCO2(ik))
        ! k in correlation is pressure based
    end do
else
    mix=1
    ! if single species is presented, by pass mixing calculation
    ! by calling the simplest sum model
endif

if (xH2O.gt.0.d0) then
    do ik=1, nq
        kp=k(ik)/xH2O
        call fskdh2oms(Tg, Tp, kp, xH2O, gH2O(ik))
        ! k in correlation is pressure based
    enddo
else
    mix=1
endif

select case (mix)
    case (1) ! only one species no mixing model is used. 
            !sum is sufficient, because either gH2O of gCO2 are all 1
        g=gH2O+gCO2-1.d0
    case (3) ! M&R mix model
        g=MRmix(k, gH2O, k, gCO2, k, 2*nq, 0) 
    case (2) ! multiplication 
        g=gH2O*gCO2
    case default
        errflag=1
        print*, "Unknown mixing model number"
end select

end subroutine fsk_corr_mix


!-----------------------------------------------------------------------------------

subroutine afun(grr,krr, grc,krc, g, w,a)
real(dp), dimension(:), intent(in) :: grr, grc, g
real(dp), dimension(size(grr,1)), intent(in) :: krr
real(dp), dimension(size(grc,1)), intent(in) :: krc
real(dp), dimension(size(g,1)), intent(in) :: w
real(dp), dimension(size(g,1)), intent(out) :: a
! internal
real(dp), dimension(0:size(g,1)) :: gb, gq, kq ! bounds
integer :: iq, nq, nrr,nrc
nq = size(g,1); nrr = size(grr,1); nrc = size(grc,1)
gb(0)=0.d0
do iq=1,nq
    gb(iq) = gb(iq-1)+w(iq)
enddo
call linearInterpMono(nrr, grr, krr, nq, gb(1:), kq(1:))
call linearInterpMono(nrc, krc, grc, nq, kq(1:), gq(1:))
gq(0)=0.d0
a=(gq(1:nq)-gq(0:nq-1))/w
end subroutine

end module libJCai


