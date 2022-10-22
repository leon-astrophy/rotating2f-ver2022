!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutines initialize the computational grid !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETGRID
USE DEFINITION
IMPLICIT NONE

! Dummy integers !
INTEGER :: i, j, k

! Assign distance !
DO j = 1, NDIV
	r(j) = rmax * (DBLE(j) - 1.0D0)/(DBLE(NDIV) - 1.0D0)
END DO
DO j = 1, KDIV
	mu(j) = (DBLE(j) - 1.0D0)/(DBLE(KDIV) - 1.0D0)
END DO

! Polar distance !
DO i = 1, KDIV
	DO j = 1, NDIV
		r_polar(i, j) = r(j)*SQRT(1.0D0 - mu(i)**2)
	END DO
END DO

END SUBROUTINE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine guesses the initial density distributions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE INITIALRHO
USE DEFINITION
IMPLICIT NONE

! Dummy integers !
INTEGER :: i, j, k

! For DM !
IF(DM_flag == 1) THEN

	! Assign uniform density at the beginning of iteration 
	DO j = 1, NDIV
		If(j <= ra1) THEN
			rho1(:,j) = -(r(j)/r(ra1))**2 + 1.0D0
		ELSE	
			rho1(:,j) = 0.0D0
		END IF
	END DO

END IF

! For NM !
DO j = 1, NDIV
	If(j <= ra2) THEN
		rho2(:,j) = -(r(j)/r(ra2))**2 + 1.0D0
	ELSE	
		rho2(:,j) = 0.0D0
	END IF
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the scaling constants for ideal fermi gas EOS !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETCONST
USE DEFINITION
IMPLICIT NONE

! Find the scaling constant !
a_max1 = (me1**4*clight**5)/(2.4D1*pi**2*h_bar**3)
b_max1 = ((mb1/ye1)*me1**3*clight**3)/(3.0D0*pi**2*h_bar**3)

! NM !
b_e = (me2**3*clight**3)/(3.0D0*pi**2*h_bar**3)
a_max2 = (me2**4*clight**5)/(2.4D1*pi**2*h_bar**3)
b_max2 = ((mb2/ye2 + me2)*me2**3*clight**3)/(3.0D0*pi**2*h_bar**3)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine gets the parameters of the parameterized EOS !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETPARAM
USE DEFINITION
IMPLICIT NONE

! Get constant !
IF(n13_flag == 1) THEN
	rho_1 = 2.0D7
	rho_2 = 2.0D13	
	y_1 = 0.5D0
	y_2 = 0.285D0
	y_c = 0.035D0
ELSEIF(g15_flag == 1) THEN
	rho_1 = 3.0D7
	rho_2 = 2.0D13	
	y_1 = 0.5D0
	y_2 = 0.278D0
	y_c = 0.035D0
ELSEIF(s15_flag == 1) THEN
	rho_1 = 2.2D8
	rho_2 = 9.5D12
	y_1 = 0.5D0
	y_2 = 0.279D0
	y_c = 0.022D0
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initializes the axial boundary point for the NM !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BOUNDARY_NM
USE DEFINITION
IMPLICIT NONE

! Dummy !
REAL (DP) :: dummy

! NM Index !
dummy = (DBLE(NDIV) - 1.0D0)/rmax + 1.0D0
ra2 = NINT(dummy)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine finds the legendre polynominal !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETLEGENDRE
USE DEFINITION 
IMPLICIT NONE

! Dummy integers !
INTEGER :: i, j, k

! Find all the legendre function !
DO i = 1, KDIV
	pn(i,0) = 1.0D0
	pn(i,1) = mu(i)
	DO j = 2, LMAX*2
		CALL LEGENDRE(pn(i,j), mu(i), j)
	END DO
END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine returns the legendre function !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LEGENDRE(pout, x, in)
USE DEFINITION 
IMPLICIT NONE

! Input integer !
INTEGER, INTENT(IN) :: in

! Input real number !
REAL (DP), INTENT(IN) :: x

! Output polynominal !
REAL (DP), INTENT(OUT) :: pout

! Dummy integer !
INTEGER :: j

! legender function !
REAL (DP), DIMENSION(0:2*LMAX) :: pnx

! Assign !
pnx(0) = 1.0D0
pnx(1) = x

! For higher legender function !
DO j = 2, in
	pnx(j) = ((2.0D0*DBLE(j) - 1.0D0)*pnx(j-1)*x - (DBLE(j) - 1.0D0)*pnx(j-2))/DBLE(j)
END DO

pout = pnx(in)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine does a mulitpole expansion !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GETEXPANSION
USE DEFINITION 
IMPLICIT NONE

! Dummy integers !
INTEGER :: i, j, k, m, l

! Do the loop !
DO m = 0, LMAX*2
	DO j = 1, NDIV
		DO k = 1, NDIV
			fp(k,j,m) = fpoisson(r(k),r(j),DBLE(m))
		END DO
	END DO
END DO

contains
	real(DP) function fpoisson(rprime,rin,index)
	implicit none
	real(DP) :: rprime,rin,index
	If (rprime<rin) THEN
		fpoisson = rprime**(index+2.0D0)/rin**(index+1.0D0)
	ELSEIF (rprime>rin) THEN
		fpoisson = rin**(index)/rprime**(index-1.0D0)
	ELSEIF (rprime == rin) THEN
		fpoisson = rprime
	END IF
	end function

END SUBROUTINE