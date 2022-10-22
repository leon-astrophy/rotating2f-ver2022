!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine computes the rotational potentials !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FINDROTATION
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! Initialize !
psi2 = 0.0D0

! Find rotational potential !
If(rigid2 == 1) THEN
	DO i = 1, KDIV
		DO j = 1, NDIV
			psi2(i,j) = rigid_potential(r_polar(i, j))
		END DO
	END DO
ELSEIF(vconst2 == 1) THEN
	DO i = 1, KDIV
		DO j = 1, NDIV
			psi2(i,j) = vconst_potential(r_polar(i, j))
		END DO
	END DO
ELSEIF(jconst2 == 1) THEN
	DO i = 1, KDIV
		DO j = 1, NDIV
			psi2(i,j) = jconst_potential(r_polar(i, j))
		END DO
	END DO
ELSEIF(kepler2 == 1) THEN
	DO i = 1, KDIV
		DO j = 1, NDIV
			psi2(i,j) = kepler_potential(r_polar(i, j))
		END DO
	END DO
ELSEIF(yoon2 == 1) THEN
	DO i = 1, KDIV
		DO j = 1, NDIV
			psi2(i,j) = rigid_potential(r_polar(i, j))
			IF(r_polar(i, j) > rcore) THEN
				psi2(i,j) = psi2(i,j) - c0*(yoon_potential (r_polar(i, j)) - yoon_potential (rcore))
			END IF
		END DO
	END DO
ELSEIF(awd2 == 1) THEN
	DO i = 1, KDIV
		DO j = 1, NDIV
			psi2(i,j) = rigid_potential(r_polar(i, j))
			IF(r_polar(i, j) > rcore) THEN
				psi2(i,j) = psi2(i,j) - awd_potential (r_polar(i, j))
			END IF
		END DO
	END DO
END IF

! Adjust the potential !
DO j = 1, NDIV
	IF(r(j) > rout) THEN
		psi2(:,j) = 0.0D0
	END IF
END DO

! For DM rotation !
IF(dm_flag == 1) THEN

	! No DM rotation !
	psi1 = 0.0D0

END IF

contains

	real(DP) function rigid_potential (r)
	implicit none
	real(DP) :: r
	rigid_potential = - 0.5D0*r*r
	end function

	real(DP) function vconst_potential (r)
	implicit none
	real(DP) :: r
	vconst_potential = - 0.5D0*log(r*r + dconst*dconst)
	end function

	real(DP) function jconst_potential (r)
	implicit none
	real(DP) :: r
	jconst_potential = 0.5D0/(r*r + dconst*dconst)
	end function

	real(DP) function kepler_potential (r)
	implicit none
	real(DP) :: r, rootd, rootr, part1, part2, part3, part4
	rootd = sqrt(dconst)
	rootr = sqrt(r)
	part1 = -6.0D0*rootr/(rootd**3 + rootr**3)
	part2 = 2.0D0*log(rootd + rootr)/dconst
	part3 = -log(dconst + r - rootd*rootr)/dconst
	part4 = -2.0D0*sqrt(3.0D0)*atan((1.0D0 - 2.0D0*rootr/rootd)/sqrt(3.0D0))/dconst
	kepler_potential = -(1.0D0/9.0D0)*(part1 + part2 + part3 + part4)
	end function

	real(DP) function yoon_potential (r)
	implicit none
	real(DP) :: r, part1, part2
	part1 = 5.0D0*atan(sqrt((r - rcore)/rcore))/(6.4D1*rcore**(1.5D0))
	part2 = sqrt(r - rcore)*(-4.8D1*rcore**3 + 1.36D2*rcore**2*r - 1.18D2*rcore*r**2 + 1.5D1*r**3)/(1.92D2*rcore*r**4)
	yoon_potential = part1 + part2
	end function

	real(DP) function awd_potential (r)
	implicit none
	real(DP) :: r, part1, part2, part3, part4, part5, part6, part7, part8, part9, part10
	part1 = 0.25D0*r**4*(c1**2 + 2.0D0*c2)
	part2 = r**6*(2.0D0*c1*c3 + c2**2 + 2.0D0*c4)/6.0D0
	part3 = 0.125D0*r**8*(2.0D0*c1*c5 + 2.0D0*c2*c4 + c3**2)
	part4 = 2.0D0*r**7*(c1*c4 + c2*c3 + c5)/7.0D0
	part5 = 0.4D0*r**5*(c1*c2 + c3)
	part6 = 2.0D0*r**3*c1/3.0D0
	part7 = 2.0D0*r**9*(c2*c5 + c3*c4)/9.0D0
	part8 = 0.1D0*r**10*(2.0D0*c3*c5 + c4**2)
	part9 = 2.0D0*r**11*c4*c5/11.0D0
	part10 = r**12*c5**2/12.0D0
	awd_potential = part1 + part2 + part3 + part4 + part5 + part6 + part7 + part8 + part9 + part10
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine computes the rotational velocities !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINDVELOCITY
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! Initialize !
velp2 = 0.0D0

! Find rotational potential !
If(rigid2 == 1) THEN
	DO i = 1, KDIV
		DO j = 1, NDIV
			omega2(i,j) = rigid_omega2 (r_polar(i, j))
			velp2(i,j) = rigid_omega2 (r_polar(i, j))*r_polar(i, j)
		END DO
	END DO
ELSEIF(vconst2 == 1) THEN
	DO i = 1, KDIV
		DO j = 1, NDIV
			omega2(i,j) = vconst_omega2 (r_polar(i, j))
			velp2(i,j) = vconst_omega2 (r_polar(i, j))*r_polar(i, j)
		END DO
	END DO
ELSEIF(jconst2 == 1) THEN
	DO i = 1, KDIV
		DO j = 1, NDIV
			omega2(i,j) = jconst_omega2 (r_polar(i, j))
			velp2(i,j) = jconst_omega2 (r_polar(i, j))*r_polar(i, j)
		END DO
	END DO
ELSEIF(kepler2 == 1) THEN
	DO i = 1, KDIV
		DO j = 1, NDIV
			omega2(i,j) = kepler_omega2 (r_polar(i, j))
			velp2(i,j) = kepler_omega2 (r_polar(i, j))*r_polar(i, j)
		END DO
	END DO
ELSEIF(yoon2 == 1) THEN
	DO i = 1, KDIV
		DO j = 1, NDIV
			IF(r_polar(i, j) < rcore) THEN
				omega2(i,j) = rigid_omega2 (r_polar(i, j))
				velp2(i,j) = rigid_omega2 (r_polar(i, j))*r_polar(i, j)
			ELSE
				omega2(i,j) = yoon_omega2 (r_polar(i, j))
				velp2(i,j) = yoon_omega2 (r_polar(i, j))*r_polar(i, j)
			END IF
		END DO
	END DO
ELSEIF(awd2 == 1) THEN
	DO i = 1, KDIV
		DO j = 1, NDIV
			IF(r_polar(i, j) < rcore) THEN
				omega2(i,j) = rigid_omega2 (r_polar(i, j))
				velp2(i,j) = rigid_omega2 (r_polar(i, j))*r_polar(i, j)
			ELSE
				omega2(i,j) = awd_omega2 (r_polar(i, j))
				velp2(i,j) = awd_omega2 (r_polar(i, j))*r_polar(i, j)
			END IF
		END DO
	END DO
END IF

! For DM !
IF(dm_flag == 1) THEN

	! No DM rotation !
	velp1 = 0.0D0
	omega1 = 0.0D0

END IF

contains

	real(DP) function rigid_omega1 (r)
	implicit none
	real(DP) :: r
	rigid_omega1 = SQRT(h02new1)
	end function

	real(DP) function rigid_omega2 (r)
	implicit none
	real(DP) :: r
	rigid_omega2 = SQRT(h02new2)
	end function

	real(DP) function vconst_omega1 (r)
	implicit none
	real(DP) :: r
	vconst_omega1 = SQRT(h02new1)/SQRT(dconst**2 + r**2)
	end function

	real(DP) function vconst_omega2 (r)
	implicit none
	real(DP) :: r
	vconst_omega2 = SQRT(h02new2)/SQRT(dconst**2 + r**2)
	end function

	real(DP) function jconst_omega1 (r)
	implicit none
	real(DP) :: r
	jconst_omega1 = SQRT(h02new1)/(dconst**2 + r**2)
	end function

	real(DP) function jconst_omega2 (r)
	implicit none
	real(DP) :: r
	jconst_omega2 = SQRT(h02new2)/(dconst**2 + r**2)
	end function

	real(DP) function kepler_omega2 (r)
	implicit none
	real(DP) :: r
	kepler_omega2 = SQRT(h02new2)/(dconst**(1.5D0) + r**(1.5D0))
	end function

	real(DP) function yoon_omega2 (r)
	implicit none
	real(DP) :: r
	yoon_omega2 = SQRT(h02new2)*SQRT(c0*(r - rcore)**(2.5D0)/r**6 + 1.0D0)
	end function

	real(DP) function awd_omega2 (r)
	implicit none
	real(DP) :: r
	awd_omega2 = SQRT(h02new2)*(1.0D0 + c1*r + c2*r**2 + c3*r**3 + c4*r**4 + c5*r**5)
	end function

END SUBROUTINE