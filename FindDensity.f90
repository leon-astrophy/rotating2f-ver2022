!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine returns the atmospheric enthalpies in CGS units !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_ATMENTHALPY
USE DEFINITION
IMPLICIT NONE

! DM Atmospheric enthalpy !
ha1 = (8.0D0*a_max1/b_max1)

! NM Atmospheric enthalpy !
IF(fermi_flag == 1) THEN
	ha2 = (8.0D0*a_max2/b_max2)
ELSE
	ha2 = 0.0D0
	hfloor2 = 10.0D0**(htable2(0))
END IF

! DM Atmospheric density !
rhoa1 = 0.0D0

! NM Atmospheric density !
rhoa2 = 0.0D0

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine returns the maximum enthalpies in CGS units !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GET_MAXENTHALPY
USE DEFINITION
IMPLICIT NONE

! NM Atmospheric enthalpy !
IF(fermi_flag == 1) THEN
	hmax2 = (8.0D0*a_max2/b_max2)*SQRT((rhomax2/b_max2)**(2.0D0/3.0D0) + 1.0D0)
ELSE
	IF(rhomax2 > 10.0D0**(rhotable2(0))) THEN
		CALL EOSRtoH(rhomax2, hmax2)
	ELSE
		hmax2 = (8.0D0*a_max2/b_max2)*(SQRT((rhomax2/b_max2)**(2.0D0/3.0D0) + 1.0D0) - 1.0D0)
	END IF
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine returns the densities of DM/NM in code units !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINDDENSITY
USE IEEE_ARITHMETIC
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! Dummy variables !
REAL (DP) :: ye, xrho, h_scale

! Scaling !
h_scale = (hmax2/hhmaxnew2)

! Fermi EOS for NM !
IF(fermi_flag == 1) THEN
	DO i = 1, KDIV
		DO j = 1, NDIV
			IF(hhatnew2(i,j) < hhata2) THEN
				rho2(i,j) = 0.0D0
			ELSEIF(hhatnew2(i,j) == hhata2) THEN
				rho2(i,j) = rhoa2
			ELSE
				rho2(i,j) = ((b_max2*hhatnew2(i,j)*h_scale/8.0D0/a_max2)**2 - 1.0D0)**3
				IF(rho2(i,j) < 0.0D0) THEN
					rho2(i,j) = 0.0D0
				END IF
				rho2(i,j) = SQRT(rho2(i,j))
				rho2(i,j) = rho2(i,j)*(b_max2/rhomax2)
			END IF
			dlfmmo2 (i,j) = (rho2(i,j)*rhomax2/b_max2)**(1.0D0/3.0D0)
		END DO	
	END DO
ELSEIF(yex_flag == 1) THEN
	DO i = 1, KDIV
		DO j = 1, NDIV
			IF(hhatnew2(i,j) >= hhatfloor2) THEN
				CALL EOSHtoX(hhatnew2(i,j), dlfmmo2 (i,j))
				ye = 1.0D0/(2.0D0 + 1.255D-2*dlfmmo2(i,j) + 1.755D-5*dlfmmo2(i,j)**2 + 1.376D-6*dlfmmo2(i,j)**3)
				rho2(i,j) = (mb2/ye + me2)*b_e*dlfmmo2 (i,j)**3/rhomax2
			ELSEIF(hhatnew2(i,j) > hhata2) THEN
				rho2(i,j) = (((b_max2*hhatnew2(i,j)*h_scale/8.0D0/a_max2) + 1.0D0)**2 - 1.0D0)**3
				IF(rho2(i,j) < 0.0D0) THEN
					rho2(i,j) = 0.0D0
				END IF
				rho2(i,j) = SQRT(rho2(i,j))
				rho2(i,j) = rho2(i,j)*(b_max2/rhomax2)
				dlfmmo2 (i,j) = (rho2(i,j)*rhomax2/b_max2)**(1.0D0/3.0D0)
		      	ELSEIF(hhatnew2(i,j) == hhata2) THEN                                                                                                                                                                                                                                                                 
				rho2(i,j) = rhoa2
				dlfmmo2 (i,j) = (rho2(i,j)*rhomax2/b_max2)**(1.0D0/3.0D0)
			ELSE
				rho2(i,j) = 0.0D0
				dlfmmo2 (i,j) = (rho2(i,j)*rhomax2/b_max2)**(1.0D0/3.0D0)
			END IF
		END DO	
	END DO
ELSEIF(yerho_flag == 1) THEN
	DO i = 1, KDIV
		DO j = 1, NDIV
			IF(hhatnew2(i,j) >= hhatfloor2) THEN
				CALL EOSHtoR(hhatnew2(i,j), rho2 (i,j))
				IF(vul_flag == 1) THEN
					CALL EOSRtoYE(rho2 (i,j), ye)
				ELSE
					xrho = max(-1.0D0, min(1.0D0, ((2.0D0*log10(rho2(i,j)) - log10(rho_2) - log10(rho_1))/(log10(rho_2) - log10(rho_1)))))
					ye = 0.5D0*(y_2 + y_1) + 0.5D0*xrho*(y_2 - y_1) + y_c*(1.0D0 - abs(xrho) + 4.0D0*abs(xrho)*(abs(xrho) - 0.5D0)*(abs(xrho) - 1.0D0))
				END IF
				dlfmmo2 (i,j) = (rho2(i,j)/b_e/(mb2/ye + me2))**(1.0D0/3.0D0)
				rho2(i,j) = rho2(i,j)/rhomax2
			ELSEIF(hhatnew2(i,j) > hhata2) THEN
				rho2(i,j) = (((b_max2*hhatnew2(i,j)*h_scale/8.0D0/a_max2) + 1.0D0)**2 - 1.0D0)**3
				IF(rho2(i,j) < 0.0D0) THEN
					rho2(i,j) = 0.0D0
				END IF
				rho2(i,j) = SQRT(rho2(i,j))
				rho2(i,j) = rho2(i,j)*(b_max2/rhomax2)
				dlfmmo2 (i,j) = (rho2(i,j)*rhomax2/b_max2)**(1.0D0/3.0D0)
		      	ELSEIF(hhatnew2(i,j) == hhata2) THEN                                                                                                                                                                                                                                                                 
				rho2(i,j) = rhoa2
				dlfmmo2 (i,j) = (rho2(i,j)*rhomax2/b_max2)**(1.0D0/3.0D0)
			ELSE
				rho2(i,j) = 0.0D0
				dlfmmo2 (i,j) = (rho2(i,j)*rhomax2/b_max2)**(1.0D0/3.0D0)
			END IF
		END DO	
	END DO
ELSEIF(hw_flag == 1) THEN
	DO i = 1, KDIV
		DO j = 1, NDIV
			IF(hhatnew2(i,j) >= hhatfloor2) THEN
				CALL EOSHtoR(hhatnew2(i,j), rho2 (i,j))
				CALL EOSHtoX(hhatnew2(i,j), dlfmmo2 (i,j))
				rho2(i,j) = rho2(i,j)/rhomax2
			ELSEIF(hhatnew2(i,j) > hhata2) THEN
				rho2(i,j) = (((b_max2*hhatnew2(i,j)*h_scale/8.0D0/a_max2) + 1.0D0)**2 - 1.0D0)**3
				IF(rho2(i,j) < 0.0D0) THEN
					rho2(i,j) = 0.0D0
				END IF
				rho2(i,j) = SQRT(rho2(i,j))
				rho2(i,j) = rho2(i,j)*(b_max2/rhomax2)
				dlfmmo2 (i,j) = (rho2(i,j)*rhomax2/b_max2)**(1.0D0/3.0D0)
		      	ELSEIF(hhatnew2(i,j) == hhata2) THEN                                                                                                                                                                                                                                                                 
				rho2(i,j) = rhoa2
				dlfmmo2 (i,j) = (rho2(i,j)*rhomax2/b_max2)**(1.0D0/3.0D0)
			ELSE
				rho2(i,j) = 0.0D0
				dlfmmo2 (i,j) = (rho2(i,j)*rhomax2/b_max2)**(1.0D0/3.0D0)
			END IF
		END DO	
	END DO
END IF

! For DM input !
IF(dm_flag == 1) THEN

	! Fermi EOS for DM !
	DO i = 1, KDIV
		DO j = 1, NDIV
			IF(hhatnew1(i,j) < hhata1) THEN
				rho1(i,j) = 0.0D0
			ELSEIF(hhatnew1(i,j) == hhata1) THEN
				rho1(i,j) = rhoa1
			ELSE
				rho1 (i,j)= ((b_max1*hhatnew1(i,j)*(hmax1/hhmaxnew1)/8.0D0/a_max1)**2 - 1.0D0)**3
				IF(rho1(i,j) < 0.0D0) THEN
					rho1(i,j) = 0.0D0
				END IF
				rho1(i,j) = SQRT(rho1(i,j))
				rho1(i,j) = rho1(i,j)*(b_max1/rhomax2)
			END IF
			dlfmmo1 (i,j) = (rho1(i,j)*rhomax2/b_max1)**(1.0D0/3.0D0)
		END DO
	END DO

END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine returns the pressure in code units !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINDPRESSURE
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! Get NM new pressure !
DO i = 1, KDIV
	DO j = 1, NDIV
		IF(rho2(i,j) >= rhoa2) THEN
			IF(dlfmmo2(i,j) > 3.0D-3) THEN
				p2 (i,j) = a_max2*large_pressure(dlfmmo2(i,j))
			ELSE
				p2 (i,j) = a_max2*small_pressure(dlfmmo2(i,j))
			END IF
		ELSE
			p2 (i,j) = 0.0D0
		END IF 
	END DO
END DO

! Scale the pressure !
p2 = (p2/pressure)

! Get DM new pressure !
IF(dm_flag == 1) THEN
	DO i = 1, KDIV
		DO j = 1, NDIV
			IF(rho1(i,j) >= rhoa1) THEN
				IF(dlfmmo1(i,j) > 3.0D-3) THEN
					p1 (i,j) = a_max1*large_pressure(dlfmmo1(i,j))
				ELSE
					p1 (i,j) = a_max1*small_pressure(dlfmmo1(i,j))
				END IF
			ELSE
				p1 (i,j) = 0.0D0
			END IF 
		END DO
	END DO

	! Scale the pressure !
	p1 = (p1/pressure)
END IF

contains
	real(DP) function large_pressure(x)
	implicit none
	real(DP) :: x
	large_pressure = x*SQRT(x**2 + 1.0D0)*(2.0D0*x**2 - 3.0D0) + 3.0D0*log(x + SQRT(x**2 + 1.0D0))
	end function

	real(DP) function small_pressure(x)
	implicit none
	real(DP) :: x
	small_pressure = 1.6D0*x**5 - (4.0D0/7.0D0)*x**7 + (1.0D0/3.0D0)*x**9 - (5.0D0/2.2D1)*x**11 & 
			+ (3.5D1/2.08D2)*x**13 - (2.1D1/1.6D2)*x**15 + (2.31D2/2.176D3)*x**17
	end function

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine returns the internal energies in code units !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FINDEPSILON
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! Get NM new pressure !
DO i = 1, KDIV
	DO j = 1, NDIV
		IF(rho2(i,j) >= rhoa2) THEN
			IF(dlfmmo2(i,j) > 3.0D-3) THEN
				eps2 (i,j) = a_max2*large_energy(dlfmmo2(i,j))
			ELSE
				eps2 (i,j) = a_max2*small_energy(dlfmmo2(i,j))
			END IF
		ELSE
			eps2 (i,j) = 0.0D0
		END IF 
	END DO
END DO

! Scale the pressure !
eps2 = (eps2/pressure)

! Get DM new pressure !
IF(dm_flag == 1) THEN
	DO i = 1, KDIV
		DO j = 1, NDIV
			IF(rho1(i,j) >= rhoa1) THEN
				IF(dlfmmo1(i,j) > 3.0D-3) THEN
					eps1 (i,j) = a_max1*large_energy(dlfmmo1(i,j))
				ELSE
					eps1 (i,j) = a_max1*small_energy(dlfmmo1(i,j))
				END IF
			ELSE
				eps1 (i,j) = 0.0D0
			END IF 
		END DO
	END DO

	! Scale the pressure !
	eps1 = (eps1/pressure)
END IF

contains
	real(DP) function large_energy(x)
	implicit none
	real(DP) :: x
	large_energy = 3.0D0*x*SQRT(x**2 + 1.0D0)*(1.0D0 + 2.0D0*x**2) - 3.0D0*log(x + SQRT(x**2 + 1.0D0)) - 8.0D0*x**3
	end function

	real(DP) function small_energy(x)
	implicit none
	real(DP) :: x
	small_energy = 8.0D0*x**3 + (1.2D1/5.0D0)*x**5 - (3.0D0/7.0D0)*x**7 + (1.0D0/6.0D0)*x**9 - (1.5D1/1.76D2)*x**11 & 
			+ (2.1D1/4.16D2)*x**13 - (2.1D1/6.40D2)*x**15 + (9.9D1/4.352D3)*x**17 - 8.0D0*x**3
	end function

END SUBROUTINE