!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This program generates rotating stars using the self-consistent, iterative method !
! c.f. Hachisu 1986 apj, 61:479-507	  					    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM ROTATING
USE DEFINITION
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE (*,*) '-----------------------'
WRITE (*,*) '|Rotating Star Program|'
WRITE (*,*) '-----------------------'
WRITE (*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Read EOS table !
IF(fermi_flag == 1) THEN
	ye2 = 5.0E-1_DP
ELSE
	CALL EOSTABLE
	IF(yerho_flag == 1) THEN
		IF(vul_flag == 1) THEN
			CALL YETABLE
		ELSE
			CALL GETPARAM
		END IF
	END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Initialize computational boxes !
CALL GETGRID

! Get fermi gas EOS constants !
CALL GETCONST

! Get atmospheric enthalpies !
CALL GET_ATMENTHALPY

! Get legendre polynominals !
CALL GETLEGENDRE

! Get poisson equations multipole expansions !
CALL GETEXPANSION

! Find rotational potential !
CALL FINDROTATION

! Get NM boundary !
CALL BOUNDARY_NM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Start computing model !
IF(dm_flag == 1) THEN
	WRITE (*,*) 'Start Solving For Two-Fluid Star'
	WRITE (*,*)
	CALL SCF2F
ELSE
	WRITE (*,*) 'Start Solving For One-Fluid Star'
	WRITE (*,*)
	CALL SCF1F
END IF

! Output parameter files for data analysis
CALL OUTPUT_PARAM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Print out !
WRITE (*,*) 'Done!'
 
END PROGRAM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine check if the model reaches the critical rotational limit !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CHECKMODEL
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k, r1_max

! Check equatorial radius !
DO j = 1, NDIV
	IF(rho2(1,j) > rhoa2 .AND. rho2(1,j+1) == rhoa2) THEN
		r_equator2 = r(j+1)
		EXIT 
	END IF
END DO

! Check axis radius !
DO j = 1, NDIV
	IF(rho2(KDIV,j) > rhoa2 .AND. rho2(KDIV,j+1) == rhoa2) THEN
		r_axis2 = r(j+1)
		EXIT 
	END IF
END DO

! Determine whether the model is successful !
IF(rb2 == 1) THEN
	IF(r_equator2 /= r(ra2)) THEN
		success_flag = .false.
	END IF
ELSE	
	IF(r_equator2 /= r(ra2) .OR. r_axis2 /= r(rb2)) THEN
		success_flag = .false.
	END IF
END IF

! For DM !
r1_max = 0
IF(DM_flag == 1) THEN
	DO j = 1, NDIV
		IF(rho1(1,j) > rhoa1 .AND. rho1(1,j+1) == rhoa1) THEN
			r1_max = j + 1
			EXIT 
		END IF
	END DO
	IF(r1_max == 0 .OR. r1_max == NDIV) THEN
		success_flag = .false.
	END IF
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine remove unwanted NM !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CHECKENTHALPY
USE DEFINITION
IMPLICIT NONE

! Integer !
INTEGER :: i, j, k

! Modify density only for the NM component !
DO i = 1, KDIV
	DO j = 1, NDIV
		IF(r(j) > rout) THEN
			hhatnew2(i,j) = 0.0D0 !hhata2
		END IF
	END DO
END DO

END SUBROUTINE