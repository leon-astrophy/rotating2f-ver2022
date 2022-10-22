#import required package#
import os
import numpy as np
import matplotlib.pyplot as plt
import math
from decimal import Decimal, ROUND_HALF_UP
# filename#
filename = []

for root, dirs, files in os.walk('../Profile/'):
    for file in files:
        if file.endswith("NM.dat"):
            filename.append(os.path.join(file))

# loop#
for o in range(0, len(filename)):
   load = np.loadtxt('../Profile/' + filename[o], max_rows=1)
   KDIV, NDIV, rmax = int(load[0]), int(load[1]), load[2]
   rho2 = np.loadtxt('../Profile/' + filename[o], skiprows=1)

   # assign grid value to x and z axis#
   polar = np.zeros(KDIV)
   for i in range(0, KDIV):
      polar[i] = math.pi / 2 - math.acos(float(i) / (float(KDIV) - 1))
   position = np.zeros(NDIV)
   for i in range(0, NDIV):
      position[i] = rmax * (float(i) / (float(NDIV) - 1))
      # mesh the grid for x and z#
      r, theta = np.meshgrid(polar, position)

   # Generate Contour Plot #
   levels=[0.001, 0.01,0.1,0.2,0.3,0.4,0.6,0.9]
   fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
   contourplot = ax.contourf(r,theta, rho2, 10, cmap='seismic')
   cax = ax.contour(r, theta, rho2, levels, colors='white', linestyles='dashed')
   plt.colorbar(contourplot)
   ax.set_thetamin(0)
   ax.set_thetamax(90)
   ax.clabel(cax, inline=True, fontsize=10)

   #ax.set_title('NM Density Contour Plot')
   plt.grid(True)
   plt.savefig(filename[o]+'.png')
   plt.clf()

   plt.plot(position, rho2[:,0])
   plt.plot(position, rho2[:,KDIV-1])
   plt.show()
   plt.clf()

# filename#
filename = []

for root, dirs, files in os.walk('../Profile/'):
    for file in files:
        if file.endswith("DM.dat"):
            filename.append(os.path.join(file))

# loop#
for o in range(0, len(filename)):
   load = np.loadtxt('../Profile/' + filename[o], max_rows=1)
   KDIV, NDIV, rmax = int(load[0]), int(load[1]), load[2]
   rho2 = np.loadtxt('../Profile/' + filename[o], skiprows=1)

   # assign grid value to x and z axis#
   polar = np.zeros(KDIV)
   for i in range(0, KDIV):
      polar[i] = math.pi / 2 - math.acos(float(i) / (float(KDIV) - 1))
   position = np.zeros(NDIV)
   for i in range(0, NDIV):
      position[i] = rmax * (float(i) / (float(NDIV) - 1))
      # mesh the grid for x and z#
      r, theta = np.meshgrid(polar, position)

   # Generate Contour Plot #
   levels=[0.001, 0.01,0.1,0.2,0.3,0.4,0.6,0.9]
   fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
   contourplot = ax.contourf(r,theta, rho2, 10, cmap='seismic')
   cax = ax.contour(r, theta, rho2, levels, colors='white', linestyles='dashed')
   plt.colorbar(contourplot)
   ax.set_thetamin(0)
   ax.set_thetamax(90)
   ax.clabel(cax, inline=True, fontsize=10)

   #ax.set_title('NM Density Contour Plot')
   plt.grid(True)
   plt.savefig(filename[o]+'.png')
   plt.clf()

1/0

1/0
#loadeddata#
parameter=np.genfromtxt('../Profile/Star_WENO_Parameter.dat')
condition=np.genfromtxt('../Profile/Star_WENO_Parameter.dat',dtype = 'str')

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

#assign steps#
dindex = (nend - nstart)/(nindex - 1)

#assign radius#
radius2 = 1

#do the loop#
for m in range (0, nindex):

   #assign#
   ninput = nstart + m*dindex
   ninput = Decimal(ninput).quantize(Decimal('.000'), ROUND_HALF_UP)

   for n in range(axstart,axend - 1,daxratio):

      #assign#
      axispt = n
      axis2 = rmax*(axispt - 1)/(NDIV - 1)
      axratio = axis2/radius2
      axratio = Decimal(axratio).quantize(Decimal('.000'), ROUND_HALF_UP)

      #loaddata#
      if (dmflag == 1):
         rho1=np.loadtxt('../Profile/Star_WENO_Density_NMCentralDensity_'+"{:.3f}".format(ninput)+'_AxisRatio_'"{:.3f}".format(axratio)+'_DMMass_'"{:.3f}".format(dmmass)+'_DM.dat',skiprows=1)
         rho2=np.loadtxt('../Profile/Star_WENO_Density_NMCentralDensity_'+"{:.3f}".format(ninput)+'_AxisRatio_'"{:.3f}".format(axratio)+'_DMMass_'"{:.3f}".format(dmmass)+'_NM.dat',skiprows=1)
      else :
         rho2=np.loadtxt('../Profile/Star_WENO_Density_NMCentralDensity_'+"{:.3f}".format(ninput)+'_AxisRatio_'"{:.3f}".format(axratio)+'_NM.dat',skiprows=1)

      #assign grid value to x and z axis#
      polar = np.zeros(KDIV)
      for i in range (0, KDIV):
         polar [i] = math.pi/2 - math.acos(float(i)/(float(KDIV) - 1))
      position = np.zeros(NDIV)
      for i in range (0, NDIV):
         position[i] = rmax*(float(i)/(float(NDIV) - 1))

      # mesh the grid for x and z#
      r, theta = np.meshgrid(polar, position)

      # Generate Contour Plot #
      if (dmflag == 1):
         levels = [0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.6, 0.9]
         fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
         contourplot = ax.contourf(r, theta, rho1/np.max(rho1), 50, cmap='viridis')
         cax = ax.contour(r, theta, rho1/np.max(rho1), levels, colors='white', linestyles='dashed')
         plt.colorbar(contourplot)
         ax.set_thetamin(0)
         ax.set_thetamax(90)
         ax.clabel(cax, inline=True, fontsize=10)
         #ax.set_title('Normalized DM Density Contour Plot')
         plt.grid(True)
         plt.savefig('DMDensityContour_NMCentralDensity_'+"{:.3f}".format(ninput)+'_AxisRatio_'"{:.3f}".format(axratio)+'_DMMass_'"{:.3f}".format(dmmass)+'.png')
         plt.clf()

      # Generate Contour Plot #
      levels=[0.001, 0.01,0.1,0.2,0.3,0.4,0.6,0.9]
      fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
      contourplot = ax.contourf(r,theta, rho2, 50, cmap='seismic')
      cax = ax.contour(r, theta, rho2, levels, colors='white', linestyles='dashed')
      plt.colorbar(contourplot)
      ax.set_thetamin(0)
      ax.set_thetamax(90)
      ax.clabel(cax, inline=True, fontsize=10)
      #ax.set_title('NM Density Contour Plot')
      plt.grid(True)
      if (dmflag == 1):
         plt.savefig('NMDensityContour_NMCentralDensity_'+"{:.3f}".format(ninput)+'_AxisRatio_'"{:.3f}".format(axratio)+'_DMMass_'"{:.3f}".format(dmmass)+'.png')
      else:
         plt.savefig('NMDensityContour_NMCentralDensity_'+"{:.3f}".format(ninput)+'_AxisRatio_'"{:.3f}".format(axratio)+'.png')
      plt.clf()

      #######################################################################################################################################

      if (dmflag == 1):

         plt.plot(position,rho1[:,0]/np.max(rho1), label='Polar')
         plt.plot(position,rho1[:,KDIV-1]/np.max(rho1), label='Equatorial',linestyle='--')
         plt.xlabel('Dimensionless Distance')
         plt.ylabel('Dimensionless Density')
         plt.title('Normalized Polar And Equatorial DM Density Profile')
         plt.grid(True)
         plt.legend(loc='upper right')
         plt.savefig('DMDensityProfile_NMCentralDensity_'+"{:.3f}".format(ninput)+'_AxisRatio_'"{:.3f}".format(axratio)+'_DMMass_'"{:.3f}".format(dmmass)+'.png')
         plt.clf()

      plt.plot(position,rho2[:,0], label='Polar')
      plt.plot(position,rho2[:,KDIV-1], label='Equatorial',linestyle='--')
      plt.xlabel('Dimensionless Distance')
      plt.ylabel('Dimensionless Density')
      plt.title('Polar And Equatorial NM Density Profile')
      plt.grid(True)
      plt.legend(loc='upper right')
      if (dmflag == 1):
         plt.savefig('NMDensityProfile_NMCentralDensity_'+"{:.3f}".format(ninput)+'_AxisRatio_'"{:.3f}".format(axratio)+'_DMMass_'"{:.3f}".format(dmmass)+'.png')
      else:
         plt.savefig('NMDensityProfile_NMCentralDensity_'+"{:.3f}".format(ninput)+'_AxisRatio_'"{:.3f}".format(axratio)+'.png')
      plt.clf()

      #prevent run-time error#
      plt.close('all')

