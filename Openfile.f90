!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine open files for densities profiles !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE OPENFILE_GLOBAL
USE DEFINITION 
IMPLICIT NONE

! Open file !
OPEN (UNIT = 21, FILE = './Parameter/Star_WENO_Global_NMCentralDensity_'//trim(str(rhoscale2))//'_NM.dat', STATUS = 'REPLACE')

! For DM !
IF(dm_flag == 1) THEN
	OPEN (UNIT = 22, FILE = './Parameter/Star_WENO_Global_NMCentralDensity_'//trim(str(rhoscale2))//'_DM.dat', STATUS = 'REPLACE')
END IF

! Header !
WRITE (21, *) '------------------------------------------------------------------------------------------------------------'
WRITE (21, *) ' masst, m2out, reout, j2out, h0out, wa2, va2, stab, bind, axratio2, ibar, qbar, rbar, rhomax2, tidal, viral '
WRITE (21, *) '------------------------------------------------------------------------------------------------------------'

! For DM !
IF(dm_flag == 1) THEN
	WRITE (22, *) '----------------------------------------------'
	WRITE (22, *) ' masst, m1out, rhoscale1, re1, stab, axratio1 '
	WRITE (22, *) '----------------------------------------------'
END IF

contains

	character(len=20) function str(k)
    	REAL (DP), INTENT(IN) :: k
    	write (str, '(f10.3)') k
    	str = adjustl(str)
	end function str

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine close files for densities profiles !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CLOSEFILE_GLOBAL
USE DEFINITION 
IMPLICIT NONE

! Close file !
CLOSE (21)

! For DM !
IF(dm_flag == 1) THEN
	CLOSE (22)
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine open files for densities profiles !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE OPENFILE_PROFILE
USE DEFINITION
IMPLICIT NONE

! Open file !
OPEN (UNIT = 41, FILE = './Profile/Star_WENO_Density_NMCentralDensity_'//trim(str(rhoscale2))//'_AxisRatio_'//trim(str(axratio2))//'_NM.dat', STATUS = 'REPLACE')

! For DM !
IF(dm_flag == 1) THEN
	OPEN (UNIT = 51, FILE = './Profile/Star_WENO_Density_NMCentralDensity_'//trim(str(rhoscale2))//'_AxisRatio_'//trim(str(axratio2))//'_DM.dat', STATUS = 'REPLACE')
END IF

contains

	character(len=20) function str(k)
    	REAL (DP), INTENT(IN) :: k
    	write (str, '(f10.3)') k
    	str = adjustl(str)
	end function str

END SUBROUTINE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine close files for densities profiles !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CLOSEFILE_PROFILE
USE DEFINITION
IMPLICIT NONE

! Close file !
CLOSE (41)

! For DM !
IF(dm_flag == 1) THEN
	CLOSE (51)
END IF

END SUBROUTINE