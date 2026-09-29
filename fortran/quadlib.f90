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
!**********************************************************************************************************
!
!**********************************************************************************************************
MODULE quadlib
CONTAINS

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
SUBROUTINE quadgen(Cheb2,g,w,NoP) !quadrature closed at 0
  IMPLICIT NONE
  LOGICAL,INTENT(IN)    :: Cheb2
  INTEGER,INTENT(IN)    :: NoP
  REAL(8),INTENT(OUT)   :: g(NoP),w(NoP)
!local variables
  REAL(8) :: gg(2*NoP-1),ww(2*NoP-1),dummy
  INTEGER(2) :: i

    CALL gaucheb2(gg,ww,2*NoP-1)
    g= -gg(NoP:2*NoP-1)
    w= ww(NoP:2*NoP-1)
    g(1)=abs(g(1)); w(1)=w(1)/2.d0
    IF (.not.Cheb2) THEN                ! transformation on Gauss-Chebyshev2 quadrature
      DO i=1, NoP
        w(i)= 1.5d0*DSQRT(1.d0-g(i))*w(i)
        g(i)= 1.d0-(1.d0-g(i))**1.5
      ENDDO
    ENDIF
    dummy= SUM(w(1:NoP))
    w(1:NoP)=w(1:NoP)/dummy

END SUBROUTINE quadgen
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
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
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
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
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

SUBROUTINE clencur(x,w,n)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx> By: Anquan Wang <xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
! DISCRIPTION:                                                                                       !
!   Calculate the zeros and weights of Clenshaw-Curtis rule                                          !
!   The zeros can be used as abscissas of Gaussian quadrature which are distributed between [-1,1]   !
!   The integral can be evaluated in the form:                                                       !
!                         Integral(f(x)*dx)= Sum(Wi*f(xi))                                           !
!----------------------------------------------------------------------------------------------------!
! PARAMETERS:                                                                                        !
!   n    -> the degree of Tchebyscheff polynomial of 2nd kind                                        !
!   x(n) -> the array of abscissas                                                                   !
!   w(n) -> the array of weights                                                                     !
!----------------------------------------------------------------------------------------------------!
! MEMO:                                                                                              !
!   P. J. Davis and P. Rabinowitz, Methods of Numerical Integration, 2nd edt, pp86                   !
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx> 2003-06-02 <xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
IMPLICIT NONE
  INTEGER :: n
  DOUBLE PRECISION,PARAMETER :: PI= 3.1415926535897932384626
  REAL(8) :: x(n),w(n)
  DOUBLE PRECISION :: theta, sum
  INTEGER :: m,k
  x(1)=1.d0
  x(n)=-1.d0
  IF (MOD(n,2).eq.1) THEN  ! J.Cai add ".eq.1" to eliminate ifort ext
    w(1)= 1.d0/REAL(n*(n-2)) ! odd
  ELSE
    W(1)= 1.d0/REAL((n-1)*(n-1)) ! even
  ENDIF
  w(n)=w(1)
  DO k= 2,(n+1)/2
    theta= REAL(k-1)*PI/(n-1.d0)
    x(k)= DCOS(theta)
    sum= 0.d0
    DO m= 1,(n-1)/2
      sum= sum + 2.d0*DCOS(2*m*theta)/(4.d0*m*m-1.d0)
    ENDDO
    IF (MOD(n,2).eq.1)THEN ! J.Cai add ".eq.1" to eliminate ifort ext
      sum= sum - DCOS(2*m*theta)/(4.d0*m*m-1.d0)
    ENDIF
    w(k)= 1.d0 - sum
    x(n+1-k)= -x(k)
    w(n+1-k)= w(k)
  ENDDO
END SUBROUTINE clencur
!*****************************************************************************
subroutine ms_quadgen(alpha,nkg,fsk,fsg,nq,gq,wq)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!NOTE:::::
! RECOMMENDED FOR MULTI-SCALE CALCUATION TO ACCOUNT FOR OVERLAP PARAETER
! this routine calculate the quadrature points and weights
! The g space (0,1) is divided into two part, one from 0 to g(0.01)
! and the other from g(0.01) to 1. 
! The reason to have two part g is to account for lambda-k relation.
! for large optical thickness, there are more points in the first part
! of the g, and for small optical thickness, there are more points
! in the second part of g.
! Alpha is the stretch factor and is according to optical thickness
! the larger the alpha, the more points at large g's
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 implicit none
 integer,intent(in) :: nkg,nq
 real(8),dimension(nkg),intent(in) :: fsk, fsg
 real(8),intent(in) :: alpha
 real(8),dimension(nq),intent(out) :: gq, wq

 integer :: i, iq, igq, nqg
 integer :: nqg1, nqg2
 real(8) :: g01  ! 2 for the two gases
 real(8),dimension(nq) :: gqs, wqs

! initialize
 gq= 0.0d0; wq= 0.0d0
 gqs=0.0d0; wqs=0.0d0

! find g values corresponding to k= 0.01
 do i=1, nkg
    if( fsk(i) >= 1.0d-02 ) then
       g01= fsg(i)
       exit
    end if
 end do

! determine number of quad point for the two g parts
! g01 portion of total nq for the first g part, no less than 5
 nqg1 = nq * g01
 nqg2 = nq - nqg1
 if ( nqg2 < 6 ) then
    nqg2 = 6
    nqg1 = nq -nqg2
 end if

! determine the quad abscissa and weight
! fisrt g part
 nqg = nqg1
 if (nqg>0) then
    ! a Gaussian quadtrature scheme 
    call quadgen2(.false., gqs(1:nqg), wqs(1:nqg), nqg, alpha)
    ! change interval: see eq. 8 of ch. 7.2 in Kincaid $ Cheney
    do iq=1, nqg
       gq(iq)=gqs(iq)*g01
       wq(iq)=wqs(iq)*g01
    end do ! iq
!!$ write(*,*) 'fisrt g part, g='
!!$ write(*,*) (gq(iq), iq=1, nqg)
 end if
! second g part
 gqs=0.0d0; wqs=0.0d0
 nqg = nqg2
 call quadgen2(.false., gqs(1:nqg), wqs(1:nqg), nqg, alpha)
 do iq=1, nqg
    igq=nqg1 + iq
    gq(igq)=gqs(iq)*(1.0d0-g01) + g01
    wq(igq)=wqs(iq)*(1.0d0-g01)
 end do ! iq
!!$ write(6,*) 'second g part, g='
!!$ write(6,*) (gq(nqg1+iq), iq=1, nqg)

 return
end subroutine ms_quadgen
!*****************************************************************************
   SUBROUTINE gauleg(x1,x2,x,w,n)
   IMPLICIT NONE
   INTEGER           :: n
   DOUBLE PRECISION  :: x1,x2,x(n),w(n)
   DOUBLE PRECISION,PARAMETER  :: EPS=3.d-14
   INTEGER           :: i,j,m
   DOUBLE PRECISION  :: p1,p2,p3,pp,xl,xm,z,z1
   m=(n+1)/2
   xm=0.5d0*(x2+x1)
   xl=0.5d0*(x2-x1)
   do 12 i=1,m
     z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1    continue
       p1=1.d0
       p2=0.d0
       do 11 j=1,n
     p3=p2
     p2=p1
     p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11     continue
       pp=n*(z*p1-p2)/(z*z-1.d0)
       z1=z
       z=z1-p1/pp
     if(abs(z-z1).gt.EPS)goto 1
     x(i)=xm-xl*z
     x(n+1-i)=xm+xl*z
     w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
     w(n+1-i)=w(i)
12 CONTINUE
   RETURN
   END subroutine gauleg
!*****************************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
RECURSIVE SUBROUTINE QsortC(A)  ! recursive implementation of quick sort algorithm
  REAL(8), INTENT(INOUT):: A(:)
!local variables
  INTEGER :: iq

  IF(SIZE(A) > 0) THEN
     CALL Partition(A, iq)
     CALL QsortC(A(:iq-1))
     CALL QsortC(A(iq+1:))
  ENDIF
END SUBROUTINE QsortC

SUBROUTINE Partition(A, marker)
  REAL(8), INTENT(in out), DIMENSION(:) :: A
  INTEGER, INTENT(out) :: marker
  INTEGER :: i, j
  REAL(8) :: temp
  REAL(8) :: x      !pivot point
  x = A(1)
  i= 0
  j= SIZE(A) + 1
  
  DO
     j = j-1
     DO
        IF (A(j) <= x) EXIT
        j = j-1
     END DO
     i = i+1
     DO
        IF (A(i) >= x) EXIT
        i = i+1
     END DO
     IF (i < j) THEN
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     ELSE
        marker = j
        RETURN
     ENDIF
  END DO
  
END SUBROUTINE Partition


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
SUBROUTINE SPLMI (N, X, Y, B, C, D) ! 1-D monotonic spline interpolation: coefficients
      INTEGER :: N
      REAL(8) :: X(N), Y(N), B(N), C(N), D(N),DELTA(N),H(N),AL(N),BE(N)
!
!  THE COEFFICIENTS B(I), C(I), AND D(I), I=1,2,...,N ARE COMPUTED
!  FOR A MONOTONICALLY VARYING CUBIC INTERPOLATING SPLINE
!
!    S(X) = Y(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2 + D(I)*(X-X(I))**3
!
!    FOR  X(I) .LE. X .LE. X(I+1)
!    WITH Y(I+1).GE.Y(I) (ALL I) OR Y(I+1).LT.Y(I) (ALL I)
!
!  INPUT..
!
!    N = THE NUMBER OF DATA POINTS OR KNOTS (N.GE.2)
!    X = THE ABSCISSAS OF THE KNOTS IN STRICTLY INCREASING ORDER
!    Y = THE ORDINATES OF THE KNOTS
!
!  OUTPUT..
!
!    B, C, D  = ARRAYS OF SPLINE COEFFICIENTS AS DEFINED ABOVE.
!
!*************************************************************************
!
! THEORY FROM 'MONOTONE PIECEWISE CUBIC INTERPOLATION',
! BY F.N. FRITSCH AND R.E. CARLSON IN SIAM J.NUMER.ANAL.,V.17,P.238
!
!*************************************************************************
!
!    Y(I) = S(X(I))
!    B(I) = SP(X(I))
!    C(I) = SPP(X(I))/2
!    D(I) = SPPP(X(I))/6  (DERIVATIVE FROM THE RIGHT)
!
!  THE ACCOMPANYING FUNCTION SUBPROGRAM  SEVAL  CAN BE USED
!  TO EVALUATE THE SPLINE.
!
!
      NM1 = N-1
      IF ( N .LT. 2 ) RETURN
      IF ( N .LT. 3 ) GO TO 100
!
! CALCULATE THE H(I) AND DELTA(I)
!
      DO 10 I=1,NM1
      AL(I)=0.
      BE(I)=0.
      H(I)=X(I+1)-X(I)
   10 DELTA(I)=(Y(I+1)-Y(I))/H(I)
!
! CALCULATE FIRST VALUES FOR AL AND BE BY 3-POINT DIFFERENCE
!
      IF(DELTA(1).EQ.0) GOTO 15
      AL(1)=((H(1)+H(2))**2*Y(2)-H(1)**2*Y(3)-H(2)*(2.*H(1)+H(2))&
             *Y(1))/(H(2)*(H(1)+H(2))*(Y(2)-Y(1)))
   15 DO 20 I=2,NM1
      IF(DELTA(I).EQ.0) GOTO 20
      AL(I)=(H(I-1)**2*Y(I+1)+(H(I)**2-H(I-1)**2)*Y(I)-H(I)**2*&
             Y(I-1))/(H(I-1)*(H(I)+H(I-1))*(Y(I+1)-Y(I)))
   20 CONTINUE
!
      NM2=N-2
      DO 30 I=1,NM2
      IF(DELTA(I).EQ.0.) GOTO 30
      BE(I)=(H(I)**2*Y(I+2)+(H(I+1)**2-H(I)**2)*Y(I+1)-H(I+1)**2*&
             Y(I))/(H(I+1)*(H(I)+H(I+1))*(Y(I+1)-Y(I)))
   30 CONTINUE
!
      IF(DELTA(N-1).EQ.0.) GOTO 35
      BE(N-1)=(H(N-2)*(2.*H(N-1)+H(N-2))*Y(N)-(H(N-1)+H(N-2))**2&
             *Y(N-1)+H(N-1)**2*Y(N-2))/(H(N-2)*(H(N-1)+H(N-2))&
             *(Y(N)-Y(N-1)))
!
! CORRECT VALUES FOR AL AND BE
!
   35 DO 40 I=1,NM1
      IF(AL(I)+BE(I).LE.2.) GOTO 40
      IF(2.*AL(I)+BE(I).LE.3.) GOTO 40
      IF(AL(I)+2.*BE(I).LE.3.) GOTO 40
      PHI=AL(I)-(2.*AL(I)+BE(I)-3.)**2/(AL(I)+BE(I)-2.)/3.
      IF(PHI.GE.0.) GOTO 40
      TI=3./SQRT(AL(I)**2+BE(I)**2)
      AL(I)=TI*AL(I)
      BE(I)=TI*BE(I)
   40 CONTINUE
!
! CALCULATE SPLINE COEFFICIENTS
!
      DO 50 I=1,NM1
      D(I)=(AL(I)+BE(I)-2.)*DELTA(I)/(H(I)**2)
      C(I)=(3.-2.*AL(I)-BE(I))*DELTA(I)/(H(I))
      B(I)=AL(I)*DELTA(I)
      IF(B(I)*DELTA(I).GE.0.) GOTO 50
        !write(*,*)b(i),delta(i)
      B(I)=DELTA(I)
      C(I)=0.
      D(I)=0.
   50 CONTINUE
      B(N)=BE(N-1)*DELTA(N-1)
      C(N)=0.
      D(N)=0.
      RETURN
!
  100 B(1) = (Y(2)-Y(1))/(X(2)-X(1))
      C(1) = 0.
      D(1) = 0.
      B(2) = B(1)
      C(2) = 0.
      D(2) = 0.
      RETURN
END SUBROUTINE SPLMI
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      SUBROUTINE SPLINE (N, X, Y, B, C, D)
      INTEGER ::  N
      REAL(8) :: X(N), Y(N), B(N), C(N), D(N)
!
!  THE COEFFICIENTS B(I), C(I), AND D(I), I=1,2,...,N ARE COMPUTED
!  FOR A CUBIC INTERPOLATING SPLINE
!
!    S(X) = Y(I) + B(I)*(X-X(I)) + C(I)*(X-X(I))**2 + D(I)*(X-X(I))**3
!
!    FOR  X(I) .LE. X .LE. X(I+1)
!
!  INPUT..
!
!    N = THE NUMBER OF DATA POINTS OR KNOTS (N.GE.2)
!    X = THE ABSCISSAS OF THE KNOTS IN STRICTLY INCREASING ORDER
!    Y = THE ORDINATES OF THE KNOTS
!
!  OUTPUT..
!
!    B, C, D  = ARRAYS OF SPLINE COEFFICIENTS AS DEFINED ABOVE.
!
!  USING  P  TO DENOTE DIFFERENTIATION,
!
!    Y(I) = S(X(I))
!    B(I) = SP(X(I))
!    C(I) = SPP(X(I))/2
!    D(I) = SPPP(X(I))/6  (DERIVATIVE FROM THE RIGHT)
!
!  THE ACCOMPANYING FUNCTION SUBPROGRAM  SEVAL  CAN BE USED
!  TO EVALUATE THE SPLINE.
!
!
      INTEGER NM1, IB, I
      REAL T
!
      NM1 = N-1
      IF ( N .LT. 2 ) RETURN
      IF ( N .LT. 3 ) GO TO 50
!
!  SET UP TRIDIAGONAL SYSTEM
!
!  B = DIAGONAL, D = OFFDIAGONAL, C = RIGHT HAND SIDE.
!
      D(1) = X(2) - X(1)
      C(2) = (Y(2) - Y(1))/D(1)
      DO 10 I = 2, NM1
         D(I) = X(I+1) - X(I)
         B(I) = 2.*(D(I-1) + D(I))
         C(I+1) = (Y(I+1) - Y(I))/D(I)
         C(I) = C(I+1) - C(I)
   10 CONTINUE
!
!  END CONDITIONS.  THIRD DERIVATIVES AT  X(1)  AND  X(N)
!  OBTAINED FROM DIVIDED DIFFERENCES
!                                                                                 
      B(1) = -D(1)
      B(N) = -D(N-1)
      C(1) = 0.
      C(N) = 0.
      IF ( N .EQ. 3 ) GO TO 15
      C(1) = C(3)/(X(4)-X(2)) - C(2)/(X(3)-X(1))
      C(N) = C(N-1)/(X(N)-X(N-2)) - C(N-2)/(X(N-1)-X(N-3))
      C(1) = C(1)*D(1)**2/(X(4)-X(1))
      C(N) = -C(N)*D(N-1)**2/(X(N)-X(N-3))
!
!  FORWARD ELIMINATION
!
   15 DO 20 I = 2, N
         T = D(I-1)/B(I-1)
         B(I) = B(I) - T*D(I-1)
         C(I) = C(I) - T*C(I-1)
   20 CONTINUE
!
!  BACK SUBSTITUTION
!
      C(N) = C(N)/B(N)
      DO 30 IB = 1, NM1
         I = N-IB
         C(I) = (C(I) - D(I)*C(I+1))/B(I)
   30 CONTINUE
!
!  C(I) IS NOW THE SIGMA(I) OF THE TEXT
!
!  COMPUTE POLYNOMIAL COEFFICIENTS
!
      B(N) = (Y(N) - Y(NM1))/D(NM1) + D(NM1)*(C(NM1) + 2.*C(N))
      DO 40 I = 1, NM1
         B(I) = (Y(I+1) - Y(I))/D(I) - D(I)*(C(I+1) + 2.*C(I))
         D(I) = (C(I+1) - C(I))/D(I)
         C(I) = 3.*C(I)
   40 CONTINUE
      C(N) = 3.*C(N)
      D(N) = D(N-1)
      RETURN
!
   50 B(1) = (Y(2)-Y(1))/(X(2)-X(1))
      C(1) = 0.
      D(1) = 0.
      B(2) = B(1)
      C(2) = 0.
      D(2) = 0.
      RETURN
      END subroutine spline
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
SUBROUTINE SEVALA(N,X,Y,B,C,D,NA,XA,YA)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: N,NA
  REAL(8),INTENT(IN) :: X(N),Y(N),B(N),C(N),D(N),XA(NA)
  REAL(8),INTENT(OUT):: YA(NA)
!!$ local variables
  INTEGER :: I,J,DIR,M,LJ,RJ,LL,LR,RL,RR
  REAL(8) :: DX
  LL= 1; LR= N; RR= N
  LJ= 1; RJ= NA
  DO WHILE(LJ<=RJ)
     IF (XA(LJ) <= X(1) ) THEN
        YA(LJ)= B(1)*(XA(LJ)-X(1))+Y(1)
     ELSEIF( XA(LJ) >= X(N)) THEN
        YA(LJ)= B(N)*(XA(LJ)-X(N))+Y(N)
     ELSE
        IF (XA(LJ)>X(LL+1)) THEN
           DO WHILE (LL+1<LR)
              M= (LL+LR)/2
              IF (XA(LJ)>X(M)) THEN
                 LL= M
              ELSE
                 LR= M
              ENDIF
           ENDDO
        ENDIF
        DX= XA(LJ)-X(LL)
        YA(LJ)= Y(LL) + DX*(B(LL) + DX*(C(LL) + DX*D(LL)))
     ENDIF
     RL= LL
                                                                                                                                                             
     IF (RJ/= LJ) THEN
        IF (XA(RJ) <= X(1) ) THEN
           YA(RJ)= B(1)*(XA(RJ)-X(1))+Y(1)
        ELSEIF( XA(RJ) >= X(N)) THEN
           YA(RJ)= B(N)*(XA(RJ)-X(N))+Y(N)
        ELSE
           IF (XA(RJ)>X(RL+1)) THEN
              DO WHILE (RL+1<RR)
                 M= (RL+RR)/2
                 IF (XA(RJ)>X(M)) THEN
                    RL= M
                 ELSE
                    RR= M
                 ENDIF
              ENDDO
           ENDIF
                                                                                                                                                             
           DX= XA(RJ)-X(RL)
           YA(RJ)= Y(RL) + DX*(B(RL) + DX*(C(RL) + DX*D(RL)))
        ENDIF
        LR= RR
     ENDIF
     LJ= LJ+1; RJ= RJ-1
  ENDDO ! WHILE
END SUBROUTINE SEVALA
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
REAL(8) FUNCTION SEVAL(N, U, X, Y, B, C, D) ! 1-D monotonic cubic spline interpolation: evaluation
      INTEGER :: N
      REAL(8) :: U, X(N), Y(N), B(N), C(N), D(N)
!
!  THIS SUBROUTINE EVALUATES THE CUBIC SPLINE FUNCTION
!
!    SEVAL = Y(I) + B(I)*(U-X(I)) + C(I)*(U-X(I))**2 + D(I)*(U-X(I))**3
!
!    WHERE  X(I) .LT. U .LT. X(I+1), USING HORNER'S RULE
!
!  IF  U .LT. X(1) THEN  I = 1  IS USED.
!  IF  U .GE. X(N) THEN  I = N  IS USED.
!
!  INPUT..
!
!    N = THE NUMBER OF DATA POINTS
!    U = THE ABSCISSA AT WHICH THE SPLINE IS TO BE EVALUATED
!    X,Y = THE ARRAYS OF DATA ABSCISSAS AND ORDINATES
!    B,C,D = ARRAYS OF SPLINE COEFFICIENTS COMPUTED BY SPLINE
!
!  IF  U  IS NOT IN THE SAME INTERVAL AS THE PREVIOUS CALL, THEN A
!  BINARY SEARCH IS PERFORMED TO DETERMINE THE PROPER INTERVAL.
!
      INTEGER I, J, K
      REAL DX
      DATA I/1/
      IF ( I .GE. N ) I = 1
      IF ( U .LT. X(I) ) GO TO 10
      IF ( U .LE. X(I+1) ) GO TO 30
!
!  BINARY SEARCH
!
   10 I = 1
      J = N+1
   20 K = (I+J)/2
      IF ( U .LT. X(K) ) J = K
      IF ( U .GE. X(K) ) I = K
      IF ( J .GT. I+1 ) GO TO 20
!
!  EVALUATE SPLINE
!
   30 DX = U - X(I)
      SEVAL = Y(I) + DX*(B(I) + DX*(C(I) + DX*D(I)))
      RETURN
      END FUNCTION

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
real(8) function sevalY(u,x,y,b,c,d)
! called by sevalint
implicit none
real(8),intent(in) :: u, x, y, b, c, d
real(8),parameter :: third=1.0d0/3.0d0
real(8) :: z,z2

z=u-x
z2=z*z
sevalY=y*z + 0.5d0*b*z2 + third*c*z2*z+0.25d0*d*z2*z2
return
end function sevalY
!-----------------------------------------------------------------------
real(8) function sevalint(n,u,x,y,b,c,d)
! perform integration using spline coefficients
implicit none
integer,intent(in) :: n
real(8),intent(in) :: u
real(8),dimension(n),intent(in) :: x, y, b, c, d

integer :: i, j, k, iend
real(8) :: bY2, bY1

sevalint=0.0d0

if ( u <= x(1) ) then
   return
else if ( u <= x(2) ) then
   iend= 1
else if ( u >= x(n) ) then
   iend=n
else
  do j=3, n
    if (u<= x(j) ) then
       i=j
       exit
    end if
  end do
  iend= i-1
end if

if (iend == 1 ) then
  k=iend
  sevalint=sevalY(u,x(k),y(k),b(k),c(k),d(k))
  
else if (iend == n ) then
  do i=1, n-1
    k=i
    by2=sevalY(x(k+1),x(k),y(k),b(k),c(k),d(k))
    sevalint= sevalint + by2
  end do
else
  do i=1, iend-1
    k=i
    bY2=sevalY(x(k+1),x(k),y(k),b(k),c(k),d(k))
    sevalint= sevalint + by2
  end do
  k=iend
  bY2=sevalY(u,x(k),y(k),b(k),c(k),d(k))
  sevalint= sevalint + by2
end if

return
end function sevalint
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      SUBROUTINE gaujac(x,w,n,alf,bet)
      INTEGER n,MAXIT
      DOUBLE PRECISION alf,bet,w(n),x(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.D-14,MAXIT=10)
      INTEGER i,its,j
      DOUBLE PRECISION alfbet,an,bn,r1,r2,r3
      DOUBLE PRECISION a,b,c,p1,p2,p3,pp,temp,z,z1
      DO 13 i=1,n
        IF(i.EQ.1)THEN
          an=alf/n
          bn=bet/n
          r1=(1.d0+alf)*(2.78d0/(4.d0+n*n)+.768d0*an/n)
          r2=1.d0+1.48d0*an+.96d0*bn+.452d0*an*an+.83d0*an*bn
          z=1.d0-r1/r2
        ELSE IF(i.EQ.2)THEN
          r1=(4.1d0+alf)/((1.d0+alf)*(1.d0+.156d0*alf))
          r2=1.d0+.06d0*(n-8.d0)*(1.d0+.12d0*alf)/n
          r3=1.d0+.012d0*bet*(1.d0+.25d0*ABS(alf))/n
          z=z-(1.d0-z)*r1*r2*r3
        ELSE IF(i.EQ.3)THEN
          r1=(1.67d0+.28d0*alf)/(1.d0+.37d0*alf)
          r2=1.d0+.22d0*(n-8.d0)/n
          r3=1.d0+8.d0*bet/((6.28d0+bet)*n*n)
          z=z-(x(1)-z)*r1*r2*r3
        ELSE IF(i.EQ.n-1)THEN
          r1=(1.d0+.235d0*bet)/(.766d0+.119d0*bet)
          r2=1.d0/(1.d0+.639d0*(n-4.d0)/(1.d0+.71d0*(n-4.d0)))
          r3=1.d0/(1.d0+20.d0*alf/((7.5d0+alf)*n*n))
          z=z+(z-x(n-3))*r1*r2*r3
        ELSE IF(i.EQ.n)THEN
          r1=(1.d0+.37d0*bet)/(1.67d0+.28d0*bet)
          r2=1.d0/(1.d0+.22d0*(n-8.d0)/n)
          r3=1.d0/(1.d0+8.d0*alf/((6.28d0+alf)*n*n))
          z=z+(z-x(n-2))*r1*r2*r3
        ELSE
          z=3.d0*x(i-1)-3.d0*x(i-2)+x(i-3)
        ENDIF
        alfbet=alf+bet
        DO 12 its=1,MAXIT
          temp=2.d0+alfbet
          p1=(alf-bet+temp*z)/2.d0
          p2=1.d0
          DO 11 j=2,n
            p3=p2
            p2=p1
            temp=2*j+alfbet
            a=2*j*(j+alfbet)*(temp-2.d0)
            b=(temp-1.d0)*(alf*alf-bet*bet+temp*(temp-2.d0)*z)
            c=2.d0*(j-1+alf)*(j-1+bet)*temp
            p1=(b*p2-c*p3)/a
11        CONTINUE
          pp=(n*(alf-bet-temp*z)*p1+2.d0*(n+alf)*(n+bet)*p2)/(temp*&
      (1.d0-z*z))
          z1=z
          z=z1-p1/pp
          IF(ABS(z-z1).LE.EPS)GOTO 1
12      CONTINUE
        !PAUSE 'too many iterations in gaujac'
		! J.Cai "PAUSE" is deleted from standard, change to a 
		! screen output with "warning"
        write (*,*) "WARNING: MAXIMUM ITERATIONS REACHED" 
1       x(i)=z
        w(i)=EXP(gammln(alf+n)+gammln(bet+n)-gammln(n+1.d0)-gammln(n+&
      alfbet+1.d0))*temp*2.d0**alfbet/(pp*p2)
13    CONTINUE
      RETURN
      END SUBROUTINE

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      FUNCTION gammln(xx)
      DOUBLE PRECISION gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,&
      24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,&
      -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*LOG(tmp)-tmp
      ser=1.000000000190015d0
      DO 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    CONTINUE
      gammln=tmp+LOG(stp*ser/x)
      RETURN
      END FUNCTION
!  (C) Copr. 1986-92 Numerical Recipes Software -)#QJ&-1Ds.      
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

END MODULE quadlib
