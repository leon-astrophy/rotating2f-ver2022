!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the Tidal Love Number !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TIDAL
USE DEFINITION
IMPLICIT NONE

! Integer !
Integer :: i, j, k, n

! Step size !
REAL (DP) :: dh

! For ODE !
REAL (DP) :: ynew, dummy

! Integrand !
REAL (DP), DIMENSION(NDIV) :: integrand

! Enclosed mass !
REAL (DP), DIMENSION(NDIV) :: mass_n

! Ratio between density !
REAL (DP), DIMENSION(NDIV) :: bigd

! Radau's equation !
REAL (DP), DIMENSION(NDIV) :: nr2

! Reducing the geometric factor
REAL (DP), DIMENSION(NDIV) :: neta2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize !
mass_n = 0.0D0
integrand = 0.0D0

! Form the integrand !
DO j = 1, NDIV
	integrand(j)  = 4.0D0*pi*r(j)**2*rho(1,j)
END DO

! Integrate Newtonian mass !
DO j = 2, NDIV
	DO i = 1, j - 1
		mass_n (j) = mass_n (j) + (1.0D0/2.0D0)*(r(i+1) - r(i))*(integrand(i+1) + integrand(i))
	END DO
END DO

! Ratio between density !
DO j = 1, NDIV
	IF(j == 1) THEN
		bigd (j) = 1.0D0
	ELSE
		bigd (j) = (4.0D0*pi*rho(1,j)*r(j)**3)/(3.0D0*mass_n (j))
	END IF
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set step size !
dh = r(2) - r(1)

! Initialize !
nr2 = 0.0D0
neta2 = 0.0D0
integrand = 0.0D0

! Integrate the ODE !
DO j = 2, NDIV
	
	! Predictor step !
	integrand(j-1) = 6.0D0*(1.0D0 - bigd (j-1)*(neta2 (j-1) + 1.0D0)) - neta2 (j-1)*(neta2 (j-1) - 2.0D0)
	nr2 (j) = nr2(j-1) + dh * integrand(j-1)
	neta2(j) = nr2(j)/r(j)
	
	! Corrector step !
	DO n = 1, 10
		dummy = 6.0D0*(1.0D0 - bigd (j)*(neta2 (j) + 1.0D0)) - neta2 (j)*(neta2 (j) - 2.0D0)
		dummy = 0.5D0*dh*(dummy + integrand(j-1))
		nr2 (j) = nr2(j-1) + dummy 
		neta2(j) = nr2(j)/r(j)
	END DO

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find DM boundary !
DO j = 1, NDIV
	IF(rho(1,j) > 0.0D0 .AND. rho(1,j+1) <= 0.0D0) THEN
		ra1 = j+1
		EXIT
	END IF
END DO

! Post process !
IF(ra1 > ra2) THEN
	ktidal = (3.0D0 - neta2(ra1))/(2.0D0*(2.0D0 + neta2(ra1)))
ELSE
	ktidal = (3.0D0 - neta2(ra2))/(2.0D0*(2.0D0 + neta2(ra2)))
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE