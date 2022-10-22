! Define double precision !
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND (15, 307)

! Total length of the computational grid !
REAL (DP), PARAMETER :: rout = 16.0D0/15.0D0
REAL (DP), PARAMETER :: rmax = 16.0D0/15.0D0

! Define grid number, they are integer multiple of 16 !
INTEGER, PARAMETER :: KDIV = 257
INTEGER, PARAMETER :: NDIV = 257

! Define the fraction of grid number to be moved !
REAL (DP), PARAMETER :: dgrid = 0.2D0

! Define the factor of decreasing step size !
REAL (DP), PARAMETER :: dratio = 0.5D0

! Error in bisections !
REAL (DP), PARAMETER :: err = 1.0D-6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Mathematical constants and physical constants !
REAL (DP), PARAMETER :: pi = 3.1415926535897932384626433832795E0_DP
REAL (DP), PARAMETER :: h_bar = 1.05457266D-27
REAL (DP), PARAMETER :: clight = 2.99792458D10
REAL (DP), PARAMETER :: gconst = 6.6743D-8
REAL (DP), PARAMETER :: solar = 1.98847D33

! Scaling constants !
REAL (DP), PARAMETER :: jscale = 1.0D50
REAL (DP), PARAMETER :: escale = 1.0D51
REAL (DP), PARAMETER :: rscale = 1.0D8

! AMU !
REAL (DP), PARAMETER :: amu = 1.66053906660D-24
REAL (DP), PARAMETER :: amumev = 931.49410242D0

! Conversion to cubic meter !
REAL (DP), PARAMETER :: cmn3tofmn3 = 1.0E-39
REAL (DP), PARAMETER :: mev2erg = 1.60217662D-6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Baryonic mass for normal matter !			
REAL (DP), PARAMETER :: mb1 = 2.0D0*1.78266191D-25

! Fermionic mass (electrons) for normal matter !			
REAL (DP), PARAMETER :: me1 = 2.0D0*1.78266191D-25

! Electron fraction for normal matter !
REAL (DP), PARAMETER :: ye1 = 1.0E0_DP

! Multiplicity factor for normal matter				
REAL (DP), PARAMETER :: gs1 = 2.0E0_DP

! Baryonic mass for normal matter !			
REAL (DP), PARAMETER :: mb2 = 1.66053906660D-24

! Fermionic mass (electrons) for normal matter !			
REAL (DP), PARAMETER :: me2 = 9.109389699D-28

! Multiplicity factor for normal matter				
REAL (DP), PARAMETER :: gs2 = 2.0E0_DP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Maximum number of legendre polynominals !
INTEGER, PARAMETER :: LMAX = 16

! Number of iteration steps !
INTEGER, PARAMETER :: nmax = 500

! Rotation law for DM !
INTEGER, PARAMETER :: rigid1 = 0
INTEGER, PARAMETER :: vconst1 = 0
INTEGER, PARAMETER :: jconst1 = 0
INTEGER, PARAMETER :: kepler1 = 0

! Rotation law for NM !
INTEGER, PARAMETER :: rigid2 = 1
INTEGER, PARAMETER :: vconst2 = 0
INTEGER, PARAMETER :: jconst2 = 0
INTEGER, PARAMETER :: kepler2 = 0
INTEGER, PARAMETER :: yoon2 = 0
INTEGER, PARAMETER :: awd2 = 0

! Constant d in rotation law !
REAL (DP), PARAMETER :: dconst = 0.1D0

! Core radius !
REAL (DP), PARAMETER :: rcore = 0.2D0

! Scaling constant for rotational profile !
REAL (DP), PARAMETER :: c0 = 100.0D0

! For awd rotation !
REAL (DP), PARAMETER :: c1 = -3.540950059D0
REAL (DP), PARAMETER :: c2 = 35.35289764D0
REAL (DP), PARAMETER :: c3 = -104.6806287D0
REAL (DP), PARAMETER :: c4 = 115.7856925D0
REAL (DP), PARAMETER :: c5 = -43.71707232D0

! Tolerance of error !
REAL (DP), PARAMETER :: tor = 1.0D-6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Number of models per maximum density !
INTEGER, PARAMETER :: n_rho = 100

! Starting log maximum density !
REAL (DP), PARAMETER :: rhostart = 6.0D0

! Ending log maximum density !
REAL (DP), PARAMETER :: rhoend = 11.1430148003D0

! Step size !
REAL (DP), PARAMETER :: drho = (rhoend - rhostart)/(DBLE(n_rho) - 1.0D0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Starting axis point !
INTEGER, PARAMETER :: axstart = NINT((DBLE(NDIV) - 1.0D0)/rmax + 1.0D0)

! Ending axis ratio !
INTEGER, PARAMETER :: axend = 1

! Number of axis point !
INTEGER, PARAMETER :: n_axis = axstart - axend + 1

! Step size !
INTEGER, PARAMETER :: daxratio = -1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Want to solve two-fluid star? !
INTEGER, PARAMETER :: dm_flag = 1

! Maximum iteration in bisection method !
INTEGER, PARAMETER :: omax = 100

! Dark matter target mass !
REAL (DP), PARAMETER :: ratio_dm = 0.05D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Skip models with critical rotations !
LOGICAL, PARAMETER :: critical_flag = .true.

! Want to output profile? !
LOGICAL, PARAMETER :: output_profile = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Use fermi gas? !
INTEGER, PARAMETER :: fermi_flag = 1

! Harrison Wheeler EOS? !
INTEGER, PARAMETER :: hw_flag = 0

! Parametrized Ye in fermi momentum EOS? !
INTEGER, PARAMETER :: yex_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Parametrized Ye in mass density EOS? !
INTEGER, PARAMETER :: yerho_flag = 0

! Which Ye in mass density EOS? !
INTEGER, PARAMETER :: n13_flag = 0
INTEGER, PARAMETER :: g15_flag = 0
INTEGER, PARAMETER :: s15_flag = 0
INTEGER, PARAMETER :: vul_flag = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!