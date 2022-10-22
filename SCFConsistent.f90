SUBROUTINE SCF1F
USE IEEE_ARITHMETIC
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Dummy integers !
INTEGER :: i, j, k, l, m, n, o

! For iterations !
REAL (DP) :: chatold2, chatnew2, fmax2

! Potentail at boundaries point !
REAL (DP) :: psia2, psib2, phia2, phib2

! Exit criteria !
REAL (DP) :: criteria1, criteria2, criteria3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the loop in generating a lot of models !
DO l = 1, n_rho
	
	! Assign NM maximum density !
	rhoscale2 = rhostart + (DBLE(l) - 1.0D0)*drho
	rhomax2 = 1.0D1**(rhoscale2)

	! Get maximum Enthalpy !
	CALL GET_MAXENTHALPY

	! openfile !
	CALL OPENFILE_GLOBAL

	! Do the loop in axis ratio !
	DO m = axstart, axend, daxratio

		! Get NM boundary !
		rb2 = m

		! Exit if the boundary is out of range !
		IF(rb2 < 1) THEN
			EXIT
		END IF

		! Find axis ratio !
		axratio2 = (r(rb2)/r(ra2))

		! Print out !
		WRITE (*,*) 'Start iterations for ...'
		WRITE (*,*) 'Central density', rhoscale2
		WRITE (*,*) 'Axis ratio', axratio2
		WRITE (*,*)
		
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Assign uniform density at the beginning of iteration 
		IF(m == axstart) THEN 
			CALL INITIALRHO
		END IF

		! Find Potential !
		CALL GRAVRHO
		CALL FINDPOTENTIAL

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Potential ar boundary point !
		phia2 = phi(1,ra2)
		phib2 = phi(KDIV,rb2)
		psia2 = psi2(1,ra2)
		psib2 = psi2(KDIV,rb2)

		! Find h02 for NM !
		IF(axratio2 == 1.0D0) THEN
			h02new2 = 0.0D0
		ELSE
			h02new2 = -(phia2 - phib2)/(psia2 - psib2)
		END IF

		! Find F = H + C !
		CALL FIND_FPOTENTIAL
		fmax2 = maxval(fhat2)		

		! Find chat !
		chatnew2 = ((ha2/hmax2)*fmax2 + phia2 + h02new2*psia2)/(1.0D0 - (ha2/hmax2))
		hhmaxnew2 = (ha2/hmax2)*(fmax2 + chatnew2)

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Initialize !
		success_flag = .true.

		! Start the loop !
		DO n = 1, nmax

			! First backup !
			h02old2 = h02new2
			chatold2 = chatnew2
			hhmaxold2 = hhmaxnew2

			! Find enthalpy !
			hhatnew2(:,:) = chatnew2 - phi(:,:) - h02new2*psi2(:,:)
			
			! Assign hhmax !
			hhmaxnew2 = maxval(hhatnew2)

			! Atmospheric enthalpy !
			hhata2 = (ha2/hmax2)*hhmaxnew2
			hhatfloor2 = (hfloor2/hmax2)*hhmaxnew2

			! Check the enthalpy !
			CALL CHECKENTHALPY

			! Get new density !
			CALL FINDDENSITY
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
			! Override enthalpy !
			hhatnew2(1,ra2) = hhata2
			hhatnew2(KDIV,rb2) = hhata2
			
			! Override density !
			rho2(1,ra2) = rhoa2
			rho2(KDIV,rb2) = rhoa2

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			! Find gravitational potential !
			CALL GRAVRHO
			CALL FINDPOTENTIAL

			! Potential ar boundary point !
			phia2 = phi(1,ra2)
			phib2 = phi(KDIV,rb2)
			
			! Find h02 !
			IF(axratio2 == 1.0D0) THEN
				h02new2 = 0.0D0
			ELSE
				h02new2 = -(phia2 - phib2)/(psia2 - psib2)
			END IF
			
			! Find modified potential !
			CALL FIND_FPOTENTIAL
			fmax2 = maxval(fhat2)

			! Find chat !
			chatnew2 = ((ha2/hmax2)*fmax2 + phia2 + h02new2*psia2)/(1.0D0 - (ha2/hmax2))

			! Exit conditions !
			IF(n > 1) THEN
				criteria1 = ABS((h02new2 - h02old2)/h02old2)
				criteria2 = ABS((chatnew2 - chatold2)/chatold2)
				criteria3 = ABS((hhmaxnew2 - hhmaxold2)/hhmaxold2)
				IF(axratio2 == 1.0D0) THEN
					IF(criteria2 < tor .AND. criteria3 < tor) EXIT
				ELSE
					IF(criteria1 < tor .AND. criteria2 < tor .AND. criteria3 < tor) EXIT
				END IF
			END IF

			! Not sucessful if the rotational velocity is negative !
			IF(h02new2 < 0.0D0) THEN
				WRITE (*,*) 'The rotational velocity is negative'
				success_flag = .false.
				EXIT 
			END IF

			! Not sucessful if the model failed to converge !
			If(n == nmax) THEN
				WRITE (*,*) 'Maximum iteration reached'
				success_flag = .false.
				EXIT
			END IF

		END DO

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Not sucessful if the model has NaN !
		IF(ieee_is_nan(h02new2) .OR. ieee_is_nan(chatnew2) .OR. ieee_is_nan(hhmaxnew2)) THEN
			WRITE (*,*) 'There is NaN'
			success_flag = .false.
		END IF

		! Be sure to check density !
		IF(success_flag .eqv. .true.) THEN
			CALL CHECKMODEL
			IF(success_flag .eqv. .false.) THEN
				WRITE (*,*) 'Mass shedding occurs'
			END IF
		END IF

		! Exit if failed to converge !
		IF(critical_flag .eqv. .true.) THEN
			IF(success_flag .eqv. .false.) THEN
				WRITE (*,*) 'Critical rotation encountered'
				WRITE (*,*) 'Move on to the next density'
				WRITE (*,*) 
				EXIT
			END IF
		END IF

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Output log !
		IF(success_flag .eqv. .true.) THEN
			CALL OUTPUT_GLOBAL
		END IF

		! Openfile !
		IF(output_profile .eqv. .true.) THEN
			CALL OPENFILE_PROFILE
			CALL OUTPUT_PROFILES
			CALL CLOSEFILE_PROFILE
		END IF

		! Print out !
		WRITE (*,*) 'Done for...'
		WRITE (*,*) 'Central density', rhoscale2
		WRITE (*,*) 'Axis ratio', axratio2
		WRITE (*,*)

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	END DO

	! Close log file !
	CALL CLOSEFILE_GLOBAL

END DO
 
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SCF2F
USE IEEE_ARITHMETIC
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Dummy integers !
INTEGER :: i, j, k, l, m, n, o

! For bisection method !
REAL (DP) :: checkm, checkmlast, ratio_temp

! For iterations !
REAL (DP) :: chatold1, chatnew1
REAL (DP) :: chatold2, chatnew2, fmax2

! Potentail at boundaries point !
REAL (DP) :: fmax2
REAL (DP) :: phic1
REAL (DP) :: psia2, psib2, phia2, phib2

! Exit criteria !
REAL (DP) :: criteria1, criteria2, criteria3
REAL (DP) :: criteria4, criteria5, criteria6

! Real !
REAL (DP) :: rho_cen, rho_old, rho_temp, drho1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Do the loop in generating a lot of models !
DO l = 1, n_rho
	
	! Assign NM maximum density !
	rhoscale2 = rhostart + (DBLE(l) - 1.0D0)*drho
	rhomax2 = 1.0D1**(rhoscale2)

	! Get maximum Enthalpy !
	CALL GET_MAXENTHALPY

	! openfile !
	CALL OPENFILE_GLOBAL

	! Do the loop in axis ratio !
	DO m = axstart, axend, daxratio

		! Get NM boundary !
		rb2 = m

		! Exit if the boundary is out of range !
		IF(rb2 < 1) THEN
			EXIT
		END IF

		! Find axis ratio !
		axratio2 = (r(rb2)/r(ra2))

		! Print out !
		WRITE (*,*) 'Start iterations for ...'
		WRITE (*,*) 'Central density', rhoscale2
		WRITE (*,*) 'Axis ratio', axratio2
		WRITE (*,*)

		! Assign !
		IF(m == axstart) THEN
			rho_cen = log10(rhomax2)
		END IF

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Bisection method
		DO o = 0, omax		

			! Assign !
			rho_temp = 10.0D0**(rho_cen)
			hmax1 = (8.0D0*a_max1/b_max1)*SQRT((rho_temp/b_max1)**(2.0D0/3.0D0) + 1.0D0)

			! Assign uniform density at the beginning of iteration
			IF(m == axstart) THEN 
				CALL INITIALRHO
			END IF

			! Find Potential !
			CALL GRAVRHO
			CALL FINDPOTENTIAL
	
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			! Potential ar boundary point !
			phic1 = phi(1,1)
			phia2 = phi(1,ra2)
			phib2 = phi(KDIV,rb2)
			psia2 = psi2(1,ra2)
			psib2 = psi2(KDIV,rb2)

			! Find h02 for NM !
			IF(axratio2 == 1.0D0) THEN
				h02new2 = 0.0D0
			ELSE
				h02new2 = -(phia2 - phib2)/(psia2 - psib2)
			END IF

			! Find F = H + C !
			CALL FIND_FPOTENTIAL
			fmax2 = maxval(fhat2)		

			! Find chat !
			chatnew2 = ((ha2/hmax2)*fmax2 + phia2 + h02new2*psia2)/(1.0D0 - (ha2/hmax2))
			hhmaxnew2 = (ha2/hmax2)*(fmax2 + chatnew2)
			chatnew1 = (hhmaxnew2/hmax2)*hmax1 + phic1
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			! Initialize !
			success_flag = .true.

			! Start the loop !
			DO n = 1, nmax

				! First backup !
				chatold1 = chatnew1
				h02old2 = h02new2
				chatold2 = chatnew2
				hhmaxold2 = hhmaxnew2

				! Find enthalpy !
				hhatnew2(:,:) = chatnew2 - phi(:,:) - h02new2*psi2(:,:)
				hhatnew1(:,:) = chatnew1 - phi(:,:)
			
				! Assign hhmax !
				hhmaxnew1 = maxval(hhatnew1)
				hhmaxnew2 = maxval(hhatnew2)

				! Atmospheric enthalpy !
				hhata1 = (ha1/hmax1)*hhmaxnew1
				hhata2 = (ha2/hmax2)*hhmaxnew2
				hhatfloor2 = (hfloor2/hmax2)*hhmaxnew2

				! Check the enthalpy !
				CALL CHECKENTHALPY

				! Get new density !
				CALL FINDDENSITY
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
				! Override enthalpy !
				hhatnew2(1,ra2) = hhata2
				hhatnew2(KDIV,rb2) = hhata2
			
				! Override density !
				rho2(1,ra2) = rhoa2
				rho2(KDIV,rb2) = rhoa2

				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				! Find gravitational potential !
				CALL GRAVRHO
				CALL FINDPOTENTIAL

				! Potential ar boundary point !
				phic1 = phi(1,1)
				phia2 = phi(1,ra2)
				phib2 = phi(KDIV,rb2)
				psia2 = psi2(1,ra2)
				psib2 = psi2(KDIV,rb2)
			
				! Find h02 !
				IF(axratio2 == 1.0D0) THEN
					h02new2 = 0.0D0
				ELSE
					h02new2 = -(phia2 - phib2)/(psia2 - psib2)
				END IF
			
				! Find modified potential !
				CALL FIND_FPOTENTIAL
				fmax2 = maxval(fhat2)

				! Find chat !
				chatnew2 = ((ha2/hmax2)*fmax2 + phia2 + h02new2*psia2)/(1.0D0 - (ha2/hmax2))
				chatnew1 = (hhmaxnew2/hmax2)*hmax1 + phic1
				
				! Exit conditions !
				IF(n > 1) THEN
					criteria1 = ABS((h02new2 - h02old2)/h02old2)
					criteria2 = ABS((chatnew2 - chatold2)/chatold2)
					criteria3 = ABS((hhmaxnew2 - hhmaxold2)/hhmaxold2)
					criteria4 = ABS((chatnew1 - chatold1)/chatold1)
					IF(axratio2 == 1.0D0) THEN
						IF(criteria2 < tor .AND. criteria3 < tor .AND. criteria4 < tor) EXIT
					ELSE
						IF(criteria1 < tor .AND. criteria2 < tor .AND. criteria3 < tor .AND. criteria4 < tor) EXIT
					END IF
				END IF
				
				! Not sucessful if the model has NaN !
				IF(ieee_is_nan(chatnew1)) THEN
					WRITE (*,*) 'There is NaN'
					success_flag = .false.
					EXIT 
				END IF

				! Not sucessful if the model has NaN !
				IF(ieee_is_nan(h02new2) .OR. ieee_is_nan(chatnew2) .OR. ieee_is_nan(hhmaxnew2)) THEN
					WRITE (*,*) 'There is NaN'
					success_flag = .false.
					EXIT 
				END IF

				! Not sucessful if the model failed to converge !	
				If(n == nmax) THEN
					WRITE (*,*) 'Maximum iteration reached'
					success_flag = .false.
					EXIT 
				END IF

			END DO

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			! Find total mass !
			CALL VOL_INT(rho1, mass1)
			CALL VOL_INT(rho2, mass2)	

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! DM mass in solar mass !
			!masstemp = mass1*(mass/solar)
			!WRITE (*,*) 'DM Mass', masstemp
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			ratio_temp = mass1/(mass1 + mass2)

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!WRITE (*,*) 'DM Ratio', ratio_temp 
			!WRITE (*,*)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Check the deviation from target mass !
			!checkm = mass_dm - masstemp
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			checkmlast = checkm
			checkm = ratio_temp - ratio_dm
			IF(o == 0) THEN
				IF(checkm > 0) THEN
					drho1 = -0.1D0
				ELSE
					drho1 = 0.1D0
				END IF
			ELSE 
				IF(n == nmax) THEN
					drho1 = -ABS(drho1)
					checkm = abs(checkm)		
				ELSE
					IF(checkmlast*checkm < 0.0D0) THEN
						drho1 = -drho1
						drho1 = drho1*0.5D0
					END IF
				END IF
			END IF
			rho_cen = rho_cen + drho1

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!WRITE (*,*) 'DM Mass', masstemp
			!WRITE (*,*) 'DM Ratio', ratio_temp
			!WRITE (*,*) 'Mass Error', ABS(checkm/mass_dm)
			!WRITE (*,*) 'Ratio Error', checkm/ratio_dm
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			! Exit !
			IF(ABS(checkm/ratio_dm) < 1.0D-4) then
				EXIT
			END IF
	
			! Convergence failure !
			IF(o == omax) THEN
				WRITE (*,*) 'Bisection limit reached'
				WRITE (*,*) 'The current DM ratio', ratio_temp
				success_flag = .false.
				EXIT 
			END IF
	
		END DO	
		
		! Backup !
		rho_old = rho_cen

		! Not sucessful if the rotational velocity is negative !
		IF(h02new2 < 0.0D0) THEN
			WRITE (*,*) 'The rotational velocity is negative'
			success_flag = .false.
		END IF

		! Be sure to check density !
		IF(success_flag .eqv. .true.) THEN
			CALL CHECKMODEL
			IF(success_flag .eqv. .false.) THEN
				WRITE (*,*) 'Mass shedding occurs'
			END IF
		END IF

		! Exit if failed to converge !
		IF(critical_flag .eqv. .true.) THEN
			IF(success_flag .eqv. .false.) THEN
				WRITE (*,*) 'Critical rotation encountered'
				WRITE (*,*) 'Move on to the next density'
				WRITE (*,*) 
				EXIT
			END IF
		END IF

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		! Output log !
		IF(success_flag .eqv. .true.) THEN
			CALL OUTPUT_GLOBAL
		END IF

		! Openfile !
		IF(output_profile .eqv. .true.) THEN
			CALL OPENFILE_PROFILE
			CALL OUTPUT_PROFILES
			CALL CLOSEFILE_PROFILE
		END IF

		! Print out !
		WRITE (*,*) 'Done for...'
		WRITE (*,*) 'Central density', rhoscale2
		WRITE (*,*) 'DM density', rho_temp
		WRITE (*,*) 'Axis ratio', axratio2
		WRITE (*,*)

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	END DO

	! Close log file !
	CALL CLOSEFILE_GLOBAL

END DO
 
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!