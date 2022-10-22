!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine output density profiles !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OUTPUT_PROFILES
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j

! For NM, Density profile !
WRITE (41, *) KDIV, NDIV, rmax, rhomax2
DO j = 1, NDIV
	WRITE (41, 701) (rho2(i,j), i = 1, KDIV)!i=KDIV, 1, -1)
END DO
WRITE (41, *)

! For DM !
IF(dm_flag == 1) THEN

	! For DM, Density  profile
	WRITE (51, *) KDIV, NDIV, rmax, rhomax2
	DO j = 1, NDIV
		WRITE (51, 701) (rho1(i,j), i = 1, KDIV)!i=KDIV, 1, -1)
	END DO
	WRITE (51, *)

END IF

! Format !
701 FORMAT (1200ES16.8)

END SUBROUTINE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine output essential quantities !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OUTPUT_GLOBAL
USE DEFINITION
IMPLICIT NONE

! Find scaling constant !
CALL GETSCALE

! Find all essential quantites
CALL FINDVELOCITY
CALL FINDPRESSURE
CALL FINDEPSILON

! Integrate essential quantities !
CALL GETQUANTITY

! Do volume integrals
CALL VOL_INT(ke2, kine2)
CALL VOL_INT(grav, gravw)
CALL VOL_INT(rho2, mass2)
CALL VOL_INT(grav2, gravw2)
CALL VOL_INT(pphi2, j2)
CALL VOL_INT(p2, pint2)
CALL VOL_INT(eps2, int2)
CALL VOL_INT(i2, itotal2)
CALL VOL_INT(q2, qtotal2)

! Do for DM !
IF(DM_flag == 1) THEN
	kine1 = 0.0D0
	CALL VOL_INT(grav1, gravw1)
	j1 = 0.0D0
	CALL VOL_INT(p1, pint1)
	CALL VOL_INT(eps1, int1)
	CALL VOL_INT(i1, itotal1)
	CALL VOL_INT(q1, qtotal1)
END IF

! Tidal love number !
CALL TIDAL

! Compute/scale the quantities !
CALL GETGLOBAL

! For relation curves !
WRITE (21, 701) masst, m2out, reout, j2out, h0out, wmax2, vmax2, stab, bind, axratio2, ibar, qbar, rbar, log10(rhomax2), ktidal, viral

! For DM !
IF(dm_flag == 1) THEN
	WRITE (22, 701) masst, m1out, rhoscale1, re1, stab, axratio1
END IF

! Format !
701 FORMAT (1200ES16.8)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine output parameter files for data analysis !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OUTPUT_PARAM
USE DEFINITION
IMPLICIT NONE

! Open !
OPEN (UNIT = 999, FILE = './Profile/Star_WENO_Parameter.dat', STATUS = 'REPLACE')

! Write !
WRITE (999,*) 'DMFlag', dm_flag
WRITE (999,*) 'DMRatio', ratio_dm
IF(rigid2 == 1) THEN
	WRITE (999,*) 'Rotation', ' rigid'
ELSEIF(vconst2 == 1) THEN
	WRITE (999,*) 'Rotation', ' vconst'
ELSEIF(jconst2 == 1) THEN
	WRITE (999,*) 'Rotation', ' jconst'
ELSEIF(kepler2 == 1) THEN
	WRITE (999,*) 'Rotation', ' kepler'
ELSEIF(yoon2 == 1) THEN
	WRITE (999,*) 'Rotation', ' yoon'
ELSEIF(awd2 == 1) THEN
	WRITE (999,*) 'Rotation', ' awd'
END IF
WRITE (999,*) 'rmax', rmax
WRITE (999,*) 'KDIV', KDIV
WRITE (999,*) 'NDIV', NDIV
WRITE (999,*) 'nrho', n_rho
WRITE (999,*) 'rhostart', rhostart
WRITE (999,*) 'rhoend', rhoend
WRITE (999,*) 'drho', drho
WRITE (999,*) 'n_axis', n_axis
WRITE (999,*) 'axstart', axstart
WRITE (999,*) 'axend', axend
WRITE (999,*) 'daxratio', daxratio

! Close !
CLOSE (999)

END SUBROUTINE