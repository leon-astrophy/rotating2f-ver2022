!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine computes the total (DM+NM) density !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GRAVRHO
USE DEFINITION
IMPLICIT NONE

! Dummy integers !
INTEGER :: i, j, k

! Add the density !
IF(DM_flag == 1) THEN
	rho (:,:) = rho1(:,:) + rho2(:,:)
ELSE
	rho (:,:) = rho2(:,:)
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine computes the gravitational potential !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINDPOTENTIAL
USE DEFINITION 
!USE OMP_LIB
IMPLICIT NONE

! Dummy integers !
INTEGER :: i, j, k, m, l

! Initialization !
d1 = 0.0D0
d2 = 0.0D0
phi = 0.0D0

! Set number of threads !
!CALL omp_set_num_threads(4)

!$OMP PARALLEL SHARED(d1, d2, phi, pn, fp, r, mu, rho, KDIV, NDIV) PRIVATE(i, j, k, m, l)

! Find d1 !
!$OMP DO COLLAPSE (2) 
DO m = 0, lmax
	DO k = 1, NDIV
		DO i = 1, KDIV-2, 2
			d1(k,m) = d1(k,m) + (1.0D0/6.0D0)*(mu(i+2) - mu(i))*(pn(i,2*m)*rho(i,k) + &
			4.0D0*pn(i+1,2*m)*rho(i+1,k) + pn(i+2,2*m)*rho(i+2,k))
		END DO
	END DO
END DO
!$OMP END DO

! Find d2 !
!$OMP DO COLLAPSE (2)
DO m = 0, lmax
	DO j = 1, NDIV
		DO k = 1, NDIV-2, 2
			d2(m,j) = d2(m,j) + (1.0D0/6.0D0)*(r(k+2) - r(k))*(fp(k,j,2*m)*d1(k,m) + &
			4.0D0*fp(k+1,j,2*m)*d1(k+1,m) + fp(k+2,j,2*m)*d1(k+2,m))
		END DO
	END DO	
END DO
!$OMP END DO

! Find potential !
!$OMP DO COLLAPSE (2)
DO i = 1, KDIV
	DO j = 1, NDIV
		DO m = 0, LMAX
			phi(i,j) = phi(i,j) + (-4.0D0*pi)*d2(m,j)*pn(i,2*m)
		END DO
	END DO
END DO
!$OMP END DO

!$OMP END PARALLEL

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine computes F = H + C !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FIND_FPOTENTIAL
USE DEFINITION
IMPLICIT NONE

! Dummy integers !
INTEGER :: i, j, k

! Initialize !
fhat2 = 0.0D0

! Assign for DM !
IF(dm_flag == 1) THEN
	fhat1 = 0.0D0
	fhat1(:,:) = - phi(:,:)
END IF

! For NM !
fhat2(:,:) = -phi(:,:) - h02new2*psi2(:,:)

END SUBROUTINE