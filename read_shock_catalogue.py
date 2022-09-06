from struct import *
import numpy as np
from collections import namedtuple 
def read_shock_list(fname):
  with open(fname, 'r') as fp:
  	data = fp.read()
  	ofs = 0
  	lng = 4
  	n   = int(unpack("i", data[ofs:ofs+lng])[0])
    #s   = 'Total number of shocks = %d' % n
    #print('Total number of shocks = ', n)

  	ofs = 4
  	lng = 8*n
  	s   = "%dl" % (n)
  	l   = np.asarray(unpack(s, data[ofs:ofs+lng]))

  	ntot = sum(l)
    #s = 'Total number of particles in shocks = %d' % ntot
    #print('Total number of particles in shocks = ', ntot)

  	ofs += lng
   	lng = 8*n
  	s   = "%dl" % (n)
  	o   = np.asarray(unpack(s, data[ofs:ofs+lng]))

   	ofs += lng
  	lng = 4*n
  	s   = "%df" % (n)
  	pd  = np.asarray(unpack(s, data[ofs:ofs+lng]))

  	ofs += lng
   	lng = 8*n
  	s   = "%dl" % (n)
  	pid = np.asarray(unpack(s, data[ofs:ofs+lng]))

  return ntot,n,l,o,pd,pid

def read_shock_data(fname,l,o):
  with open(fname,'r') as fp:
    data = fp.read()

    ntot = sum(l)

    ofs = 0
    lng = 4
    n   = int(unpack("i", data[ofs:ofs+lng])[0])
    ofs += lng
    #print(n)

    d  = np.zeros(ntot,dtype=float)
    x  = np.zeros(ntot,dtype=float)
    y  = np.zeros(ntot,dtype=float)
    z  = np.zeros(ntot,dtype=float)
    vx = np.zeros(ntot,dtype=float)
    vy = np.zeros(ntot,dtype=float)
    vz = np.zeros(ntot,dtype=float)
    ids = np.zeros(ntot,dtype=long)

    for i in range(n):

      #read in densities
      lng = 4*l[i]
      s   = "%df" % l[i]
      din = np.asarray(unpack(s, data[ofs:ofs+lng]))
      d[o[i]:o[i]+l[i]] = din[:]
      ofs += lng

      #read in x
      lng = 4*l[i]
      s   = "%df" % l[i]
      xin = np.asarray(unpack(s, data[ofs:ofs+lng]))
      x[o[i]:o[i]+l[i]] = xin[:]
      ofs += lng

      #read in y
      lng = 4*l[i]
      s   = "%df" % l[i]
      yin = np.asarray(unpack(s, data[ofs:ofs+lng]))
      y[o[i]:o[i]+l[i]] = yin[:]
      ofs += lng

      #read in z
      lng = 4*l[i]
      s   = "%df" % l[i]
      zin = np.asarray(unpack(s, data[ofs:ofs+lng]))
      z[o[i]:o[i]+l[i]] = zin[:]
      ofs += lng

      #read in vx
      lng = 4*l[i]
      s   = "%df" % l[i]
      vxin = np.asarray(unpack(s, data[ofs:ofs+lng]))
      vx[o[i]:o[i]+l[i]] = vxin[:]
      ofs += lng

      #read in vy
      lng = 4*l[i]
      s   = "%df" % l[i]
      vyin = np.asarray(unpack(s, data[ofs:ofs+lng]))
      vy[o[i]:o[i]+l[i]] = vyin[:]
      ofs += lng

      #read in vz
      lng = 4*l[i]
      s   = "%df" % l[i]
      vzin = np.asarray(unpack(s, data[ofs:ofs+lng]))
      vz[o[i]:o[i]+l[i]] = vzin[:]
      ofs += lng

      #read in ids
      lng = 8*l[i]
      s   = "%dl" % l[i]
      idsin = np.asarray(unpack(s, data[ofs:ofs+lng]))
      ids[o[i]:o[i]+l[i]] = idsin[:]
      ofs += lng

  print(d.min(),d.max())
  return n, d, x, y, z, vx, vy, vz, ids

def read_shock_data_pot(fname,l,o):
  with open(fname,'r') as fp:
    data = fp.read()

    ntot = sum(l)

    ofs = 0
    lng = 4
    n   = int(unpack("i", data[ofs:ofs+lng])[0])
    ofs += lng
    #print(n)

    d  = np.zeros(ntot,dtype=float)
    x  = np.zeros(ntot,dtype=float)
    y  = np.zeros(ntot,dtype=float)
    z  = np.zeros(ntot,dtype=float)
    vx = np.zeros(ntot,dtype=float)
    vy = np.zeros(ntot,dtype=float)
    vz = np.zeros(ntot,dtype=float)
    ids = np.zeros(ntot,dtype=long)
    pot = np.zeros(ntot,dtype=float)

    #densities
    lng = 4*ntot
    s   = "%df" % ntot
    d = np.asarray(unpack(s, data[ofs:ofs+lng]))
    ofs += lng
    #x
    lng = 4*ntot
    s   = "%df" % ntot
    x = np.asarray(unpack(s, data[ofs:ofs+lng]))
    ofs += lng
    #y
    lng = 4*ntot
    s   = "%df" % ntot
    y = np.asarray(unpack(s, data[ofs:ofs+lng]))
    ofs += lng
    #z
    lng = 4*ntot
    s   = "%df" % ntot
    z = np.asarray(unpack(s, data[ofs:ofs+lng]))
    ofs += lng
    #vx
    lng = 4*ntot
    s   = "%df" % ntot
    vx = np.asarray(unpack(s, data[ofs:ofs+lng]))
    ofs += lng
    #vy
    lng = 4*ntot
    s   = "%df" % ntot
    vy = np.asarray(unpack(s, data[ofs:ofs+lng]))
    ofs += lng
    #v
    lng = 4*ntot
    s   = "%df" % ntot
    vz = np.asarray(unpack(s, data[ofs:ofs+lng]))
    ofs += lng
    #read in ids
    lng = 8*ntot
    s   = "%dl" % ntot
    ids = np.asarray(unpack(s, data[ofs:ofs+lng]))
    ofs += lng
    #pot
    lng = 4*ntot
    s   = "%df" % ntot
    pot = np.asarray(unpack(s, data[ofs:ofs+lng]))
    ofs += lng

  return n, d, x, y, z, vx, vy, vz, ids, pot
