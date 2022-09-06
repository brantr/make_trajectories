#!/usr/bin/env python
#load some libraries
import sys
import os
from array_io import *
import numpy as np
from read_shock_catalogue import *
from calculate_shock_inertia import *

#define some numbers
snap   = 4 #snapshot

#load in the shock catalog at this snapshot time
fname_list = "peaks/peak.%04d.list.gpi" % snap
#fname_list = "peaks/peak.%04d.list" % snap
ntot,n,l,o,pd,pid = read_shock_list(fname_list)

#ntot is total number of particles at this time in dense peaks
#n is the number of peaks
#l is the length of each peak
#o is the offset of each peak
#pd is the max density of the peak
#pid is the particle id of the densest particle in the peak

fname_dat  = "peaks/peak.%04d.dat.gpi"  % snap
#fname_dat  = "peaks/peak.%04d.dat"  % snap
nin, din, xin, yin, zin, vxin, vyin, vzin, idsin = read_shock_data(fname_dat,l,o)
lin = l
oin = o


#nin is the number of particles, should equal ntot
#din are the densities
#xin, yin, zin are the positions
#vxin, vyin, vzin are the velocities
#idsin are the particle ids




#for ishock in range(0,2,1):
print n
for ishock in range(n):
#for ishock in range(0,14,1):

  #string versions
  sishock = "%08d" % ishock
  exts = "%04d." % snap




  #now we get the length and offset of the shock of interest
  #get ishock
  print("ishock = ",ishock)

  l = lin[ishock]
  o = oin[ishock]
  print("l o ids[o]",l,o,idsin[o])
  if(l>1):

    #and get the particles from this shock
    n = l
    d = din[o:o+l]
    x = xin[o:o+l]
    y = yin[o:o+l]
    z = zin[o:o+l]
    vx = vxin[o:o+l]
    vy = vyin[o:o+l]
    vz = vzin[o:o+l]
    ids = idsin[o:o+l]
  
    print("Density / x / y / z / vx / vy / vz")
    print(d.min(),d.max())
    print(x.min(),x.max())
    print(y.min(),y.max())
    print(z.min(),z.max())
    print(vx.min(),vx.max())
    print(vy.min(),vy.max())
    print(vz.min(),vz.max())

    #now we have to pick a center
    #given that we are using gaussian process
    #interpolation, we can likely identify a reasonable
    #peak.

    #let's find the densest particles that are near the peak

    #get an index to the peak
    peak_index = d.argmax()
    x_d_max = x[peak_index]
    y_d_max = y[peak_index]
    z_d_max = z[peak_index]
    #print peak_index, d[peak_index], d.max()

    #first, let's figure out about how big the shock is
    #so let's get an RMS size
    r_rms = np.mean( (x - x[peak_index])**2 + (y - y[peak_index])**2 + (z - z[peak_index])**2 )**0.5
    #print r_rms

    #define a fractional density to select particles
    n_dense_min     = 10   # set a minimum number of particles to use to define the peak location
    if(n_dense_min>l):     # possible that there are very few particles, if so -- include all
      n_dense_min = l
    dense_fraction  = 0.9  # we'll take the top decile in density
    radius_fraction = 0.25 # we'll take dense particles within r_rms/4 of the peak density

    #iteratively select dense particles until we have enough
    n_dense = 0
    while(n_dense<n_dense_min):
      dense_index = np.where( (d>dense_fraction*d.max())&(((x - x[peak_index])**2 + (y - y[peak_index])**2 + (z - z[peak_index])**2 )<(radius_fraction*r_rms)**2))[0]
      n_dense = len(dense_index) #number of dense particles used to define the center
      if((n_dense<n_dense_min)&(n_dense<l)):
        #we don't have enough dense particles
        #so let's increase radius_fraction
        #and decrease density fraction
        #and try again
        radius_fraction *= 1.1
        dense_fraction  /= 1.1

    print n_dense, n_dense_min

    #define the density-weighted center of the peak
    x_center = np.sum(x[dense_index]*d[dense_index])/np.sum(d[dense_index])
    y_center = np.sum(y[dense_index]*d[dense_index])/np.sum(d[dense_index])
    z_center = np.sum(z[dense_index]*d[dense_index])/np.sum(d[dense_index])

    #define the density-weighted velocity of the peak
    vx_center = np.sum(vx[dense_index]*d[dense_index])/np.sum(d[dense_index])
    vy_center = np.sum(vy[dense_index]*d[dense_index])/np.sum(d[dense_index])
    vz_center = np.sum(vz[dense_index]*d[dense_index])/np.sum(d[dense_index])

    #define positions and velocities relative to the peak center
    xr = x - x_center
    yr = y - y_center
    zr = z - z_center
    vxr = vx
    vyr = vy
    vzr = vz

    #now, we may want to define a new subset of
    #particles for computing the moment of inertia
    #tensor.

    #define a fractional density to select particles
    n_ixx_min     = 50   # set a minimum number of particles to use to define I_xx
    if(n_ixx_min>l):      #if there are too few particles, include all of them
      n_ixx_min = l
    dense_fraction  = 0.5 # we'll take the top half in density
    radius_fraction = 0.5 # we'll take dense particles within r_rms/2 of the peak center

    #iteratively select dense particles until we have enough
    n_ixx = 0
    print "Minimum number of particles used in principal axis calculation: ",n_ixx_min
    while(n_ixx<n_ixx_min):
      ixx_index = np.where( (d>dense_fraction*d.max())&((xr**2 + yr**2 + zr**2 )<(radius_fraction*r_rms)**2))[0]
      n_ixx = len(ixx_index) #number of dense particles used to define the center
      if((n_ixx<n_ixx_min)&(n_ixx<l)):
        #we don't have enough dense particles
        #so let's increase radius_fraction
        #and try again
        radius_fraction *= 1.1
        dense_fraction  /= 1.1
        
    print "Number of particles used in principal axis calculation: ",n_ixx

    #compute the moment of inertia tensor
    #and diagonalize it
    dw = np.zeros(l)
    dw[:] = 1.
    U,W,V = calculate_shock_inertia(d[ixx_index],xr[ixx_index],yr[ixx_index],zr[ixx_index])

    #define a rotation based on the principal axes
    axt = U[0,0]
    bxt = U[1,0]
    cxt = U[2,0]
    ayt = U[0,1]
    byt = U[1,1]
    cyt = U[2,1]
    azt = U[0,2]
    bzt = U[1,2]
    czt = U[2,2]


    #get a rotated set of positions
    x_rot = axt*xr + bxt*yr + cxt*zr
    y_rot = ayt*xr + byt*yr + cyt*zr
    z_rot = azt*xr + bzt*yr + czt*zr

    vx_center_rot = axt*vx_center + bxt*vy_center + cxt*vz_center
    vy_center_rot = ayt*vx_center + byt*vy_center + cyt*vz_center
    vz_center_rot = azt*vx_center + bzt*vy_center + czt*vz_center
    vv_center_rot = (vx_center_rot**2 + vy_center_rot**2 + vz_center_rot**2)**0.5 #magnitude of velocity

    vx_rot = axt*vxr + bxt*vyr + cxt*vzr
    vy_rot = ayt*vxr + byt*vyr + cyt*vzr
    vz_rot = azt*vxr + bzt*vyr + czt*vzr
    vv_rot = (vx_rot**2 + vy_rot**2 + vz_rot**2)**0.5


    #define shock as moving in positive x direction
    if(vx_center_rot<0):
      vx_center_rot *= -1
      vx_rot *= -1
      x_rot  *= -1


    fname = "trajectories/trajectory.single.%04d.%s.txt" % (snap,sishock)
    print("Writing ",fname)
    fp = open(fname,"w")
    s = "% 10.9e\t% 10.9e\t% 10.9e\n" % (x_center,y_center,z_center)
    fp.write(s)
    s = "% 10.9e\t% 10.9e\t% 10.9e\t% 10.9e\t% 10.9e\t% 10.9e\t% 10.9e\n" % (x_d_max,y_d_max,z_d_max,d.max(),vx_center,vy_center,vz_center)
    fp.write(s)
    s = "% 10.9e\t% 10.9e\t% 10.9e\n" % (axt,ayt,azt)
    fp.write(s)
    s = "% 10.9e\t% 10.9e\t% 10.9e\n" % (bxt,byt,bzt)
    fp.write(s)
    s = "% 10.9e\t% 10.9e\t% 10.9e\n" % (cxt,cyt,czt)
    fp.write(s)
    print(x_rot.min(),x_rot.max())
    s = "% 10.9e\t% 10.9e\n" % (5*x_rot.min(),5*x_rot.max())
    fp.write(s)
    fp.close()
