#import required package#
import numpy as np
import matplotlib.pyplot as plt
import math
from decimal import Decimal, ROUND_HALF_UP

#linear interpolation#
def linear(xin, x0, x1, y0, y1):
  output = y0 + (xin - x0)*(y1 - y0)/(x1 - x0)
  return output

#loadeddata#
parameter=np.genfromtxt('../Profile/Star_WENO_Parameter.dat')
condition=np.genfromtxt('../Profile/Star_WENO_Parameter.dat',dtype = 'str')

#speed of light #
clight = 2.99792458*10**(10)

#load parameter#
dmflag = int(parameter[0,1])
dmmass = float(parameter[1,1])
rotation = str(condition[2,1])
rmax = float(parameter[3,1])
KDIV = int(parameter[4,1])
NDIV = int(parameter[5,1])
nindex = int(parameter[6,1])
nstart = float(parameter[7,1])
nend = float(parameter[8,1])
dindex = float(parameter[9,1])
naxis = int(parameter[10,1])
axstart = int(parameter[11,1])
axend = int(parameter[12,1])
daxratio = int(parameter[13,1])

#steps#
dindex = (nend - nstart)/(nindex - 1)

#arrays#
norotate = np.ndarray(shape=(nindex,2), dtype=float)
critical = np.ndarray(shape=(nindex,2), dtype=float)

# do the loop#
for m in range(0, nindex):
    ninput = nstart + m * dindex
    ninput = Decimal(ninput).quantize(Decimal('.000'), ROUND_HALF_UP)
    #loaddata#
    if dmflag == 1 :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput)+'_DMMass_'+"{:.3f}".format(dmmass)+'_NM.dat',skiprows=3)
    else :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput) + '_NM.dat',skiprows=3)
    #array#
    masst = outputnm[:,0]
    bind = outputnm[:,8]
    size = np.size(masst)
    #non rotate curve#
    norotate[m, 0] = masst[0]
    norotate[m, 1] = bind[0]
    critical[m, 0] = masst[size - 1]
    critical[m, 1] = bind[size - 1]
    plt.plot(masst, bind, label='$\\rho_{Max}=$' + "{:.2f}".format(ninput))
l1, = plt.plot(norotate[:, 0], norotate[:, 1], label='Non-Rotation', linestyle='-.')
l2, = plt.plot(critical[:, 0], critical[:, 1], label='Critical-Rotation', linestyle='--')
plt.legend(handles=[l1, l2], loc="upper right")
plt.xlabel('Mass ($M_{\\odot}$)')
plt.ylabel('Binding Energy ($10^{51}$ erg)')
plt.title('Bindding Energy Against $M$')
plt.grid(True)
plt.savefig('bmtcurve')
#plt.show()
plt.clf()

# do the loop#
for m in range(0, nindex):
    ninput = nstart + m * dindex
    ninput = Decimal(ninput).quantize(Decimal('.000'), ROUND_HALF_UP)
    #loaddata#
    if dmflag == 1 :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput)+'_DMMass_'+"{:.3f}".format(dmmass)+'_NM.dat',skiprows=3)
    else :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput) + '_NM.dat',skiprows=3)
    #array#
    masst = outputnm[:,0]
    omega = outputnm[:,4]
    size = np.size(masst)
    #non rotate curve#
    norotate[m, 0] = masst[0]
    norotate[m, 1] = omega[0]
    critical[m, 0] = masst[size - 1]
    critical[m, 1] = omega[size - 1]
    plt.plot(masst, omega, label='$\\rho_{Max}=$' + "{:.2f}".format(ninput))
l1, = plt.plot(norotate[:, 0], norotate[:, 1], label='Non-Rotation', linestyle='-.')
l2, = plt.plot(critical[:, 0], critical[:, 1], label='Critical-Rotation', linestyle='--')
plt.legend(handles=[l1, l2], loc="upper right")
plt.xlabel('Mass ($M_{\\odot}$)')
if(rotation == 'rigid'):
   plt.ylabel('Angular Velocity (s$^{-1}$)')
   plt.title('$\\omega$ Against $M$')
if(rotation == 'vconst'):
   plt.ylabel('Velocity (cms$^{-1}$)')
   plt.title('$v$ Against $M$')
if(rotation == 'jconst'):
   plt.ylabel('Angular Momentum (cm$^{2}$s$^{-1}$)')
   plt.title('$j$ Against $M$')
if (rotation == 'kepler'):
   plt.ylabel('Kepler Rotation (cm$^{\\frac{3}{2}}$s$^{-1}$)')
   plt.title('$k$ Against $M$')
plt.grid(True)
plt.savefig('h0mtcurve')
#plt.show()
plt.clf()

# do the loop#
for m in range(0, nindex):
    ninput = nstart + m * dindex
    ninput = Decimal(ninput).quantize(Decimal('.000'), ROUND_HALF_UP)
    #loaddata#
    if dmflag == 1 :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput)+'_DMMass_'+"{:.3f}".format(dmmass)+'_NM.dat',skiprows=3)
    else :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput) + '_NM.dat',skiprows=3)
    #array#
    masst = outputnm[:,0]
    omegae = outputnm[:,5]
    size = np.size(masst)
    #non rotate curve#
    norotate[m, 0] = masst[0]
    norotate[m, 1] = omegae[0]
    critical[m, 0] = masst[size - 1]
    critical[m, 1] = omegae[size - 1]
    plt.plot(masst, omegae, label='$\\rho_{Max}=$' + "{:.2f}".format(ninput))
l1, = plt.plot(norotate[:, 0], norotate[:, 1], label='Non-Rotation', linestyle='-.')
l2, = plt.plot(critical[:, 0], critical[:, 1], label='Critical-Rotation', linestyle='--')
plt.legend(handles=[l1, l2], loc="upper right")
plt.xlabel('Mass ($M_{\\odot}$)')
plt.ylabel('Angular Velocity (rad s$^{-1}$)')
plt.title('Maximum Angular Velocity Against $M$')
plt.grid(True)
plt.savefig('omegamtcurve')
#plt.show()
plt.clf()

# do the loop#
for m in range(0, nindex):
    ninput = nstart + m * dindex
    ninput = Decimal(ninput).quantize(Decimal('.000'), ROUND_HALF_UP)
    #loaddata#
    if dmflag == 1 :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput)+'_DMMass_'+"{:.3f}".format(dmmass)+'_NM.dat',skiprows=3)
    else :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput) + '_NM.dat',skiprows=3)
    #array#
    masst = outputnm[:,0]
    vele = outputnm[:,6]
    size = np.size(masst)
    #non rotate curve#
    norotate[m, 0] = masst[0]
    norotate[m, 1] = vele[0]
    critical[m, 0] = masst[size - 1]
    critical[m, 1] = vele[size - 1]
    plt.plot(masst, vele/clight, label='$\\rho_{Max}=$' + "{:.2f}".format(ninput))
l1, = plt.plot(norotate[:, 0], norotate[:, 1]/clight, label='Non-Rotation', linestyle='-.')
l2, = plt.plot(critical[:, 0], critical[:, 1]/clight, label='Critical-Rotation', linestyle='--')
plt.legend(handles=[l1, l2], loc="upper right")
plt.xlabel('Mass ($M_{\\odot}$)')
plt.ylabel('Rotational Velocity (Speed Of Light)')
plt.title('Maximum Rotational Velocity Against $M$')
plt.grid(True)
plt.savefig('velmtcurve')
#plt.show()
plt.clf()

# do the loop#
for m in range(0, nindex):
    ninput = nstart + m * dindex
    ninput = Decimal(ninput).quantize(Decimal('.000'), ROUND_HALF_UP)
    #loaddata#
    if dmflag == 1 :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput)+'_DMMass_'+"{:.3f}".format(dmmass)+'_NM.dat',skiprows=3)
    else :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput) + '_NM.dat',skiprows=3)
    #array#
    masst = outputnm[:,0]
    reout = outputnm[:,2]
    #non rotate curve#
    norotate[m, 0] = masst[0]
    norotate[m, 1] = reout[0]
    size = np.size(masst)
    #non rotate curve#
    norotate[m, 0] = masst[0]
    norotate[m, 1] = reout[0]
    critical[m, 0] = masst[size - 1]
    critical[m, 1] = reout[size - 1]
    plt.plot(reout, masst, label='$\\rho_{Max}=$' + "{:.2f}".format(ninput))
l1, = plt.plot(norotate[:, 1], norotate[:, 0], label='Non-Rotation', linestyle='-.')
l2, = plt.plot(critical[:, 1], critical[:, 0], label='Critical-Rotation', linestyle='--')
plt.legend(handles=[l1, l2], loc="upper right")
plt.ylabel('Mass ($M_{\\odot}$)')
plt.xlabel('Equatorial Radius ($10^{8}$ cm)')
plt.title('$R_{e}$ Against $M$')
plt.grid(True)
plt.savefig('r2mtcurve')
#plt.show()
plt.clf()

# do the loop#
for m in range(0, nindex):
    ninput = nstart + m * dindex
    ninput = Decimal(ninput).quantize(Decimal('.000'), ROUND_HALF_UP)
    #loaddata#
    if dmflag == 1 :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput)+'_DMMass_'+"{:.3f}".format(dmmass)+'_NM.dat',skiprows=3)
    else :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput) + '_NM.dat',skiprows=3)
    #array#
    masst = outputnm[:,0]
    reout = outputnm[:,2]
    stab = outputnm[:, 7]
    #non rotate curve#
    norotate[m, 0] = masst[0]
    norotate[m, 1] = reout[0]
    #check stab#
    size = np.size(stab)
    #initalize#
    k = 0
    for i in range(0, size):
        if stab[i] > 0.14:
            k = i
            break
    if k > 0:
        mtmax = linear(0.14, stab[k-1], stab[k], masst[k-1], masst[k])
        remax = linear(0.14, stab[k-1], stab[k], reout[k-1], reout[k])
        critical[m, 0] = mtmax
        critical[m, 1] = remax
        mt = np.zeros(k+1)
        re2 = np.zeros(k+1)
        for i in range (0, k+1):
            if i != k:
                mt[i] = masst[i]
                re2[i] = reout[i]
            else:
                mt[i] = mtmax
                re2[i] = remax
    else:
        mt = np.zeros(size)
        re2 = np.zeros(size)
        for i in range (0, size):
            mt[i] = masst[i]
            re2[i] = reout[i]
        critical[m, 0] = mt[size - 1]
        critical[m, 1] = re2[size - 1]
    #alter#
    for i in range (0, size):
        if mt[i] == np.max(mt):
            critical[m, 0] = mt[i]
            critical[m, 1] = re2[i]
            break
    plt.plot(re2, mt, label='$\\rho_{Max}=$' + "{:.2f}".format(ninput))
l1, = plt.plot(norotate[:, 1], norotate[:, 0], label='Non-Rotation', linestyle='-.')
l2, = plt.plot(critical[:, 1], critical[:, 0], label='Critical-Rotation', linestyle='--')
plt.legend(handles=[l1, l2], loc="upper right")
plt.ylabel('Mass ($M_{\\odot}$)')
plt.xlabel('Equatorial Radius ($10^{8}$ cm)')
plt.title('$R_{e}$ Against $M$')
plt.grid(True)
plt.savefig('r2mtsecular')
#plt.show()
plt.clf()

# do the loop#
for m in range(0, nindex):
    ninput = nstart + m * dindex
    ninput = Decimal(ninput).quantize(Decimal('.000'), ROUND_HALF_UP)
    #loaddata#
    if dmflag == 1 :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput)+'_DMMass_'+"{:.3f}".format(dmmass)+'_NM.dat',skiprows=3)
    else :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput) + '_NM.dat',skiprows=3)
    #array#
    masst = outputnm[:,0]
    rhoc = outputnm[:,13]
    size = np.size(masst)
    #non rotate curve#
    norotate[m, 0] = masst[0]
    norotate[m, 1] = rhoc[0]
    critical[m, 0] = masst[size - 1]
    critical[m, 1] = rhoc[size - 1]
    plt.plot(rhoc, masst, label='$\\rho_{Max}=$' + "{:.2f}".format(ninput))
l1, = plt.plot(norotate[:, 1], norotate[:, 0], label='Non-Rotation', linestyle='-.')
l2, = plt.plot(critical[:, 1], critical[:, 0], label='Critical-Rotation', linestyle='--')
plt.legend(handles=[l1, l2], loc="lower right")
plt.ylabel('Mass ($M_{\\odot}$)')
plt.xlabel('Log10 Density (gcm$^{-3}$)')
plt.title('$\\rho_{\\max}$ Against $M$')
plt.grid(True)
plt.savefig('rho2mtcurve')
#plt.show()
plt.clf()

# do the loop#
for m in range(0, nindex):
    ninput = nstart + m * dindex
    ninput = Decimal(ninput).quantize(Decimal('.000'), ROUND_HALF_UP)
    #loaddata#
    if dmflag == 1 :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput)+'_DMMass_'+"{:.3f}".format(dmmass)+'_NM.dat',skiprows=3)
    else :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput) + '_NM.dat',skiprows=3)
    #array#
    masst = outputnm[:,0]
    rhoc = outputnm[:,13]
    stab = outputnm[:, 7]
    #non rotate curve#
    norotate[m, 0] = masst[0]
    norotate[m, 1] = rhoc[0]
    #check stab#
    size = np.size(stab)
    #initalize#
    k = 0
    for i in range(0, size):
        if stab[i] > 0.14:
            k = i
            break
    if k > 0:
        mtmax = linear(0.14, stab[k-1], stab[k], masst[k-1], masst[k])
        rhomax = linear(0.14, stab[k-1], stab[k], rhoc[k-1], rhoc[k])
        critical[m, 0] = mtmax
        critical[m, 1] = rhomax
        mt = np.zeros(k+1)
        rho2 = np.zeros(k+1)
        for i in range (0, k+1):
            if i != k:
                mt[i] = masst[i]
                rho2[i] = rhoc[i]
            else:
                mt[i] = mtmax
                rho2[i] = rhomax
    else:
        mt = np.zeros(size)
        rho2 = np.zeros(size)
        for i in range (0, size):
            mt[i] = masst[i]
            rho2[i] = rhoc[i]
        critical[m, 0] = mt[size - 1]
        critical[m, 1] = rho2[size - 1]
    #alter#
    for i in range (0, size):
        if mt[i] == np.max(mt):
            critical[m, 0] = mt[i]
            critical[m, 1] = rho2[i]
            break
    plt.plot(rho2, mt, label='$\\rho_{Max}=$' + "{:.2f}".format(ninput))
l1, = plt.plot(norotate[:, 1], norotate[:, 0], label='Non-Rotation', linestyle='-.')
l2, = plt.plot(critical[:, 1], critical[:, 0], label='Critical-Rotation', linestyle='--')
plt.legend(handles=[l1, l2], loc="lower right")
plt.ylabel('Mass ($M_{\\odot}$)')
plt.xlabel('Log10 Density (gcm$^{-3}$)')
plt.title('$\\rho_{\\max}$ Against $M$')
plt.grid(True)
plt.savefig('rho2mtsecular')
#plt.show()
plt.clf()

# do the loop#
for m in range(0, nindex):
    ninput = nstart + m * dindex
    ninput = Decimal(ninput).quantize(Decimal('.000'), ROUND_HALF_UP)
    #loaddata#
    if dmflag == 1 :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput)+'_DMMass_'+"{:.3f}".format(dmmass)+'_NM.dat',skiprows=3)
    else :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput) + '_NM.dat',skiprows=3)
    #array#
    mass2 = outputnm[:,1]
    reout = outputnm[:,2]
    size = np.size(mass2)
    #non rotate curve#
    norotate[m, 0] = mass2[0]
    norotate[m, 1] = reout[0]
    critical[m, 0] = mass2[size - 1]
    critical[m, 1] = reout[size - 1]
    plt.plot(reout, mass2, label='$\\rho_{Max}=$' + "{:.2f}".format(ninput))
l1, = plt.plot(norotate[:, 1], norotate[:, 0], label='Non-Rotation', linestyle='-.')
l2, = plt.plot(critical[:, 1], critical[:, 0], label='Critical-Rotation', linestyle='--')
plt.legend(handles=[l1, l2], loc="upper right")
plt.ylabel('Mass ($M_{\\odot}$)')
plt.xlabel('Equatorial Radius ($10^{8}$ cm)')
plt.title('$R_{e}$ Against $M_{NM}$')
plt.grid(True)
plt.savefig('r2m2curve')
#plt.show()
plt.clf()

# do the loop#
for m in range(0, nindex):
    ninput = nstart + m * dindex
    ninput = Decimal(ninput).quantize(Decimal('.000'), ROUND_HALF_UP)
    #loaddata#
    if dmflag == 1 :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput)+'_DMMass_'+"{:.3f}".format(dmmass)+'_NM.dat',skiprows=3)
    else :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput) + '_NM.dat',skiprows=3)
    #array#
    masst = outputnm[:,0]
    j2out = outputnm[:,3]
    size = np.size(masst)
    #non rotate curve#
    norotate[m, 0] = masst[0]
    norotate[m, 1] = j2out[0]
    critical[m, 0] = masst[size - 1]
    critical[m, 1] = j2out[size - 1]
    plt.plot(masst, j2out, label='$\\rho_{Max}=$' + "{:.2f}".format(ninput))
l1,=plt.plot(norotate[:, 0], norotate[:, 1], label='Non-Rotation',linestyle='-.')
l2,=plt.plot(critical[:, 0], critical[:, 1], label='Critical-Rotation',linestyle='--')
plt.legend(handles=[l1, l2], loc="upper left")
plt.xlabel('Mass ($M_{\\odot}$)')
plt.ylabel('Angular Momentum ($10^{51}$ ergs)')
plt.title('$J$ Against $M$')
plt.grid(True)
plt.savefig('j2mtcurve')
#plt.show()
plt.clf()

# do the loop#
for m in range(0, nindex):
    ninput = nstart + m * dindex
    ninput = Decimal(ninput).quantize(Decimal('.000'), ROUND_HALF_UP)
    #loaddata#
    if dmflag == 1 :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput)+'_DMMass_'+"{:.3f}".format(dmmass)+'_NM.dat',skiprows=3)
    else :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput) + '_NM.dat',skiprows=3)
    #array#
    axratio = outputnm[:,9]
    rhomax2 = outputnm[:,13]
    size = np.size(axratio)
    #non rotate curve#
    norotate[m, 0] = rhomax2[0]
    norotate[m, 1] = axratio[0]
    critical[m, 0] = rhomax2[size - 1]
    critical[m, 1] = axratio[size - 1]
    plt.plot(rhomax2, axratio, label='$\\rho_{Max}=$' + "{:.2f}".format(ninput))
l1,=plt.plot(norotate[:, 0], norotate[:, 1], label='Non-Rotation',linestyle='-.')
l2,=plt.plot(critical[:, 0], critical[:, 1], label='Critical-Rotation',linestyle='--')
plt.legend(handles=[l1, l2], loc="upper left")
plt.xlabel('Log10 Density (gcm$^{-3}$)')
plt.ylabel('Axis-Ratio')
plt.title('Axis-Ratio Against $\\rho_{2max}$')
plt.grid(True)
plt.savefig('axrho2curve')
#plt.show()
plt.clf()

# do the loop#
for m in range(0, nindex):
    ninput = nstart + m * dindex
    ninput = Decimal(ninput).quantize(Decimal('.000'), ROUND_HALF_UP)
    #loaddata#
    if dmflag == 1 :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput)+'_DMMass_'+"{:.3f}".format(dmmass)+'_NM.dat',skiprows=3)
    else :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput) + '_NM.dat',skiprows=3)
    #array#
    masst = outputnm[:,0]
    stab = outputnm[:, 7]
    size = np.size(masst)
    #non rotate curve#
    norotate[m, 0] = masst[0]
    norotate[m, 1] = stab[0]
    critical[m, 0] = masst[size - 1]
    critical[m, 1] = stab[size - 1]
    plt.plot(masst, stab, label='$\\rho_{Max}=$' + "{:.2f}".format(ninput))
if(np.amax(critical[:, 1])) > 0.14 :
    start = np.amin(norotate[:, 0])
    end = np.amax(critical[:, 0])
    secularx1 = np.linspace(start, end, num=100)
    seculary1 = np.linspace(0.14, 0.14, num=100)
    l3, = plt.plot(secularx1, seculary1, label='Secular-Limit', linestyle=':')
if(np.amax(critical[:, 0])) > 2.6 :
    start = np.amin(norotate[:, 1])
    end = np.amax(critical[:, 1])
    seculary2 = np.linspace(start, end, num=100)
    secularx2 = np.linspace(2.6, 2.6, num=100)
    l4, = plt.plot(secularx2, seculary2, label='Secular-Limit', linestyle=':')
l1,=plt.plot(norotate[:, 0], norotate[:, 1], label='Non-Rotation',linestyle='-.')
l2,=plt.plot(critical[:, 0], critical[:, 1], label='Critical-Rotation',linestyle='--')
plt.legend(handles=[l1, l2], loc="upper left")
plt.xlabel('Mass ($M_{\\odot}$)')
plt.ylabel('Stability Parameter')
plt.title('$\\frac{T}{|W|}$ Against $M$')
plt.grid(True)
plt.savefig('twmtcurve')
plt.xlim(2.5,2.7)
plt.ylim(0.13,0.15)
plt.savefig('twmtzoom')
plt.clf()

# do the loop#
for m in range(0, nindex):
    ninput = nstart + m * dindex
    ninput = Decimal(ninput).quantize(Decimal('.000'), ROUND_HALF_UP)
    #loaddata#
    if dmflag == 1 :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput)+'_DMMass_'+"{:.3f}".format(dmmass)+'_NM.dat',skiprows=3)
    else :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput) + '_NM.dat',skiprows=3)
    #array#
    itotal2 = outputnm[:,10]
    qtotal2 = outputnm[:,11]
    size = np.size(itotal2)
    #non rotate curve#
    norotate[m, 0] = qtotal2[0]
    norotate[m, 1] = itotal2[0]
    critical[m, 0] = qtotal2[size - 1]
    critical[m, 1] = itotal2[size - 1]
    plt.plot(qtotal2, itotal2, label='$\\rho_{Max}=$' + "{:.2f}".format(ninput))
l1,=plt.plot(norotate[:, 0], norotate[:, 1], label='Non-Rotation',linestyle='-.')
l2,=plt.plot(critical[:, 0], critical[:, 1], label='Critical-Rotation',linestyle='--')
plt.legend(handles=[l1, l2], loc="lower left")
plt.xlabel('Moment Of Inertia')
plt.ylabel('Mass Multipole Moment')
plt.title('Normalized $Q$ against normalized $I$')
plt.grid(True)
plt.savefig('QIcurve')
#plt.show()
plt.clf()

# do the loop#
for m in range(0, nindex):
    ninput = nstart + m * dindex
    ninput = Decimal(ninput).quantize(Decimal('.000'), ROUND_HALF_UP)
    #loaddata#
    if dmflag == 1 :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput)+'_DMMass_'+"{:.3f}".format(dmmass)+'_NM.dat',skiprows=3)
    else :
        outputnm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_'+"{:.3f}".format(ninput) + '_NM.dat',skiprows=3)
    #array#
    itotal2 = outputnm[:,10]
    qtotal2 = outputnm[:,11]
    stab = outputnm[:,7]
    #non rotate curve#
    norotate[m, 0] = qtotal2[0]
    norotate[m, 1] = itotal2[0]
    #check stab#
    size = np.size(stab)
    #initalize#
    k = 0
    for i in range(0, size):
        if stab[i] > 0.14:
            k = i
            break
    if k > 0:
        qtmax = linear(0.14, stab[k-1], stab[k], qtotal2[k-1], qtotal2[k])
        itmax = linear(0.14, stab[k-1], stab[k], itotal2[k-1], itotal2[k])
        critical[m, 0] = qtmax
        critical[m, 1] = itmax
        qt = np.zeros(k+1)
        it = np.zeros(k+1)
        for i in range (0, k+1):
            if i != k:
                qt[i] = qtotal2[i]
                it[i] = itotal2[i]
            else:
                qt[i] = qtmax
                it[i] = itmax
    else:
        qt = np.zeros(size)
        it = np.zeros(size)
        for i in range (0, size):
            qt[i] = qtotal2[i]
            it[i] = itotal2[i]
        critical[m, 0] = qt[size - 1]
        critical[m, 1] = it[size - 1]
    plt.plot(qt, it, label='$\\rho_{Max}=$' + "{:.2f}".format(ninput))
l1, = plt.plot(norotate[:, 0], norotate[:, 1], label='Non-Rotation', linestyle='-.')
l2, = plt.plot(critical[:, 0], critical[:, 1], label='Critical-Rotation', linestyle='--')
plt.legend(handles=[l1, l2], loc="lower left")
plt.xlabel('Moment Of Inertia')
plt.ylabel('Mass Multipole Moment')
plt.title('Normalized $Q$ against normalized $I$')
plt.grid(True)
plt.savefig('QIsecular')
# plt.show()
plt.clf()

#for dm#
if dmflag == 1 :

    # do the loop#
    for m in range(0, nindex):
        ninput = nstart + m * dindex
        ninput = Decimal(ninput).quantize(Decimal('.000'), ROUND_HALF_UP)
        # loaddata#
        outputdm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_' + "{:.3f}".format(ninput) + '_DMMass_' + "{:.3f}".format(dmmass) + '_DM.dat', skiprows=3)
        # array#
        masst = outputdm[:, 0]
        rhom1 = outputdm[:, 2]
        size = np.size(masst)
        # non rotate curve#
        norotate[m, 0] = masst[0]
        norotate[m, 1] = rhom1[0]
        critical[m, 0] = masst[size - 1]
        critical[m, 1] = rhom1[size - 1]
        plt.plot(masst, rhom1, label='$\\rho_{Max}=$' + "{:.2f}".format(ninput))
    l1, = plt.plot(norotate[:, 0], norotate[:, 1], label='Non-Rotation', linestyle='-.')
    l2, = plt.plot(critical[:, 0], critical[:, 1], label='Critical-Rotation', linestyle='--')
    plt.legend(handles=[l1, l2], loc="upper right")
    plt.xlabel('Mass ($M_{\\odot}$)')
    plt.ylabel('Log 10 DM Maximum Density (gcm$^{-3}$)')
    plt.title('DM Maximum Density Against $M$')
    plt.grid(True)
    plt.savefig('rho1mtcurve')
    #plt.show()
    plt.clf()

    # do the loop#
    for m in range(0, nindex):
        ninput = nstart + m * dindex
        ninput = Decimal(ninput).quantize(Decimal('.000'), ROUND_HALF_UP)
        # loaddata#
        outputdm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_' + "{:.3f}".format(ninput) + '_DMMass_' + "{:.3f}".format(dmmass) + '_DM.dat', skiprows=3)
        # array#
        masst = outputdm[:, 0]
        rad1 = outputdm[:, 3]
        size = np.size(masst)
        # non rotate curve#
        norotate[m, 0] = masst[0]
        norotate[m, 1] = rad1[0]
        critical[m, 0] = masst[size - 1]
        critical[m, 1] = rad1[size - 1]
        plt.plot(rad1, masst, label='$\\rho_{Max}=$' + "{:.2f}".format(ninput))
    l1, = plt.plot(norotate[:, 1], norotate[:, 0], label='Non-Rotation', linestyle='-.')
    l2, = plt.plot(critical[:, 1], critical[:, 0], label='Critical-Rotation', linestyle='--')
    plt.legend(handles=[l1, l2], loc="upper right")
    plt.ylabel('Mass ($M_{\\odot}$)')
    plt.xlabel('DM Equatorial Radius ($10^{8}$ cm)')
    plt.title('DM Equatorial Radius Against $M$')
    plt.grid(True)
    plt.savefig('r1mtcurve')
    #plt.show()
    plt.clf()

    # do the loop#
    for m in range(0, nindex):
        ninput = nstart + m * dindex
        ninput = Decimal(ninput).quantize(Decimal('.000'), ROUND_HALF_UP)
        # loaddata#
        outputdm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_' + "{:.3f}".format(ninput) + '_DMMass_' + "{:.3f}".format(dmmass) + '_DM.dat', skiprows=3)
        # array#
        masst = outputdm[:, 0]
        axratio1 = outputdm[:, 5]
        size = np.size(masst)
        # non rotate curve#
        norotate[m, 0] = masst[0]
        norotate[m, 1] = axratio1[0]
        critical[m, 0] = masst[size - 1]
        critical[m, 1] = axratio1[size - 1]
        plt.plot(masst, axratio1, label='$\\rho_{Max}=$' + "{:.2f}".format(ninput))
    l1, = plt.plot(norotate[:, 0], norotate[:, 1], label='Non-Rotation', linestyle='-.')
    l2, = plt.plot(critical[:, 0], critical[:, 1], label='Critical-Rotation', linestyle='--')
    plt.legend(handles=[l1, l2], loc="upper right")
    plt.xlabel('Mass ($M_{\\odot}$)')
    plt.ylabel('Axis Ratio')
    plt.title('DM Axis Ratio Against $M$')
    plt.grid(True)
    plt.savefig('ax1mtcurve')
    #plt.show()
    plt.clf()

    # do the loop#
    for m in range(0, nindex):
        ninput = nstart + m * dindex
        ninput = Decimal(ninput).quantize(Decimal('.000'), ROUND_HALF_UP)
        # loaddata#
        outputdm = np.loadtxt('../Parameter/Star_WENO_Global_NMCentralDensity_' + "{:.3f}".format(ninput) + '_DMMass_' + "{:.3f}".format(dmmass) + '_DM.dat', skiprows=3)
        # array#
        masst = outputdm[:, 0]
        mass1 = (outputdm[:, 1] - dmmass)/dmmass
        #plot#
        plt.plot(masst, mass1, label='$\\rho_{Max}=$' + "{:.2f}".format(ninput))
    plt.xlabel('Total Mass ($M_{\\odot}$)')
    plt.ylabel('Percentage Error Of DM Mass ($M_{\\odot}$)')
    plt.title('Consistency Check')
    plt.grid(True)
    plt.savefig('m1mtcurve')
    #plt.show()
    plt.clf()

#prevent run-time error#
plt.close('all')