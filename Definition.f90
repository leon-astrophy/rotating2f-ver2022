MODULE DEFINITION
IMPLICIT NONE
INCLUDE "Parameter.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Define direction arrays !
REAL (DP), DIMENSION(NDIV) :: r
REAL (DP), DIMENSION(KDIV) :: mu

! Polar distance !
REAL (DP), DIMENSION(KDIV,NDIV) :: r_polar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Gravitational potential !
REAL (DP), DIMENSION(KDIV,NDIV) :: phi
REAL (DP), DIMENSION(KDIV,NDIV) :: phi1
REAL (DP), DIMENSION(KDIV,NDIV) :: phi2

! Rotational potential !
REAL (DP), DIMENSION(KDIV,NDIV) :: psi1
REAL (DP), DIMENSION(KDIV,NDIV) :: psi2

! Density !
REAL (DP), DIMENSION(KDIV,NDIV) :: rho
REAL (DP), DIMENSION(KDIV,NDIV) :: rho1
REAL (DP), DIMENSION(KDIV,NDIV) :: rho2

! Pressure !
REAL (DP), DIMENSION(KDIV,NDIV) :: p1
REAL (DP), DIMENSION(KDIV,NDIV) :: p2

! Energy per unit volume !
REAL (DP), DIMENSION(KDIV,NDIV) :: eps1
REAL (DP), DIMENSION(KDIV,NDIV) :: eps2

! Rotational velocity !
REAL (DP), DIMENSION(KDIV,NDIV) :: velp1
REAL (DP), DIMENSION(KDIV,NDIV) :: velp2

! Angular velocity !
REAL (DP), DIMENSION(KDIV,NDIV) :: omega1
REAL (DP), DIMENSION(KDIV,NDIV) :: omega2

! Kinetic energy !
REAL (DP), DIMENSION(KDIV,NDIV) :: ke1
REAL (DP), DIMENSION(KDIV,NDIV) :: ke2

! Gravitational energy !
REAL (DP), DIMENSION(KDIV,NDIV) :: grav
REAL (DP), DIMENSION(KDIV,NDIV) :: grav1
REAL (DP), DIMENSION(KDIV,NDIV) :: grav2

! Angular momentum !
REAL (DP), DIMENSION(KDIV,NDIV) :: pphi1
REAL (DP), DIMENSION(KDIV,NDIV) :: pphi2

! Dimensonless Fermi Momentum !
REAL (DP), DIMENSION(KDIV,NDIV) :: dlfmmo1
REAL (DP), DIMENSION(KDIV,NDIV) :: dlfmmo2

! Enthalpy !
REAL (DP), DIMENSION(KDIV,NDIV) :: hhatnew1
REAL (DP), DIMENSION(KDIV,NDIV) :: hhatnew2

! Integration constant !
REAL (DP), DIMENSION(KDIV,NDIV) :: fhat1
REAL (DP), DIMENSION(KDIV,NDIV) :: fhat2

! Moment of inertia !
REAL (DP), DIMENSION(KDIV,NDIV) :: i1
REAL (DP), DIMENSION(KDIV,NDIV) :: i2

! Mass quadrapole moment !
REAL (DP), DIMENSION(KDIV,NDIV) :: q1
REAL (DP), DIMENSION(KDIV,NDIV) :: q2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Intermediate arrays for multipole expansion !
REAL (DP), DIMENSION(NDIV,0:LMAX) :: d1
REAL (DP), DIMENSION(0:LMAX,NDIV) :: d2

! Legendre function !
REAL (DP), DIMENSION(KDIV,0:LMAX*2) :: pn

! The 1/r expansion factor !
REAL (DP), DIMENSION(NDIV,NDIV,0:LMAX*2) :: fp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Number of lines in EOS table !
INTEGER :: nlines1, nlines2, yelines2

! For DM Table !
REAL (DP), ALLOCATABLE, DIMENSION(:) :: htable1
REAL (DP), ALLOCATABLE, DIMENSION(:) :: xtable1
REAL (DP), ALLOCATABLE, DIMENSION(:) :: rhotable1

! For NM Table !
REAL (DP), ALLOCATABLE, DIMENSION(:) :: htable2
REAL (DP), ALLOCATABLE, DIMENSION(:) :: xtable2
REAL (DP), ALLOCATABLE, DIMENSION(:) :: yetable2
REAL (DP), ALLOCATABLE, DIMENSION(:) :: rhotable2
REAL (DP), ALLOCATABLE, DIMENSION(:) :: rhoyetable2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Data analysis !
REAL (DP) :: gravw, stab, kine, ktidal
REAL (DP) :: gravw2, j2, j2out, kine2, pint2, int2, wmax2, vmax2
REAL (DP) :: mass2, m2out, volume2, itotal2, qtotal2
REAL (DP) :: gravw1, j1, j1out, kine1, pint1, int1, wmax1, vmax1, re1
REAL (DP) :: mass1, m1out, volume1, itotal1, qtotal1
REAL (DP) :: masst, bind, reout, h0out, rbar, qbar, ibar, jtotal, vtotal, viral

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Electron fractions !
REAL (DP) :: ye2

! Multiplicative constants !
REAL (DP) :: b_e

! For Fermi momentum !
REAL (DP) :: a_max1, b_max1
REAL (DP) :: a_max2, b_max2

! For enthalpy !
REAL (DP) :: hmax1, hhmaxold1, hhmaxnew1
REAL (DP) :: hmax2, hhmaxold2, hhmaxnew2

! Dimensionless atmospheric value !
REAL (DP) :: ha1, ha2
REAL (DP) :: rhoa1, rhoa2
REAL (DP) :: hhata1, hhata2
REAL (DP) :: hfloor1, hfloor2
REAL (DP) :: hhatfloor1, hhatfloor2

! Scaling constant !
REAL (DP) :: hconst1, hconst2

! Angular velocity !
REAL (DP) :: h02old1, h02new1
REAL (DP) :: h02old2, h02new2

! Equatorial and polar radius !
REAL (DP) :: r_equator1, r_axis1
REAL (DP) :: r_equator2, r_axis2

! Array index at boundary point !
INTEGER :: ra1, rb1
INTEGER :: ra2, rb2
 
! Maximum Density !
REAL (DP) :: rhomax1, rhoscale1
REAL (DP) :: rhomax2, rhoscale2

! Some scaling constant !
REAL (DP) :: roverb1, roverb2
REAL (DP) :: hasquare1, hasqaure2
REAL (DP) :: rbonethird1, rbonethird2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Target DM mass !
REAL (DP) :: targetdm

! Size ratio !
INTEGER :: dsize
REAL (DP) :: sizeratio

! Absolute error !
REAL (DP) :: abserror

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Equatorial radius !
REAL (DP) :: r_equatorial

! Axis ratio !
REAL (DP) :: axratio1
REAL (DP) :: axratio2

! Scaling constants !
REAL (DP) :: vol
REAL (DP) :: mass
REAL (DP) :: length
REAL (DP) :: energy
REAL (DP) :: pressure
REAL (DP) :: momentum
REAL (DP) :: potential
REAL (DP) :: clighthat
REAL (DP) :: persecond

! For sucessful model !
LOGICAL :: success_flag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For parameterized electron fractions !
REAL (DP) :: rho_1
REAL (DP) :: rho_2
REAL (DP) :: y_1
REAL (DP) :: y_2
REAL (DP) :: y_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE