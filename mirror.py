#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 22:10:20 2017

@author: winokur
"""
import numpy as np
from mw_library import tilt_calc
from mw_library import read_txyz
from mw_library import write_xyz
from mw_library import write_txyz
#sign =  1.0
sign = -1.0
tfile = '/home/winokur/dtu_projects/NaT2mw.txyz'
l2_num, aax,aay,aaz,atype, astring, abt, uc, header = read_txyz(tfile)
# Tidy up the location because Tinker's default is a bit ugly
# From the old structure to the new
aax *= sign
aaz *= sign
#uc[4]*= sign
xyz = np.array([aax,aay,aaz],dtype=float)
xyz = np.transpose(xyz)
#   Recenter chains outside unit cell periphery
#toler = np.array([0.03,0.08,0.05])
#xyz = olig_recenter(mol_seq,uc,xyz,toler)
#write_txyz_2('wtest.xyz_2',atype,astring,xyz,abt,uc,header)    
#
write_txyz('mirror.txyz',l2_num,atype,astring,aax,aay,aaz,abt,uc,header)
a =[sign*uc[0],0.,0.]
b =[0.,uc[1],0.]
c =[(uc[2]*np.cos(np.radians(uc[4]))),0.,(uc[2]*np.sin(np.radians(uc[4])))]
T2 = -0.5*np.asarray([sign*(uc[0]+uc[2]*np.cos(np.radians(uc[4]))),uc[1],sign*uc[2]*np.sin(np.radians(uc[4]))])
bs = np.asarray(b)
T2 = T2 + b #- 1.*np.asarray(c)

#xy=0.; yz=0.; xz=0.
xy,yz,xz = tilt_calc(a,b,c)
Lx=0.5*uc[0]; Ly=0.5*uc[1]; Lz=0.5*uc[2]
# Lx=25.; Ly=25.; Lz=25. # Oversized to eliminate VdW across grain boundaries
wrapx=np.zeros(l2_num); wrapy=np.zeros(l2_num); wrapz=np.zeros(l2_num)
xyz_new=np.zeros([l2_num,3]); xyz_wp=np.zeros([l2_num,3],dtype=int)
for i in range(l2_num):
    zmin=-Lz
    zmax= Lz
#    print i,xyz[i]
#    print i, aax[i],aay[i],aaz[i],xmin,xmax,ymin,ymax,zmin,zmax 
    for j in range(3):  # Maximum number of wraps expected with monoclinic symmetry
        xyz_new[i,0] = xyz[i,0]+wrapx[i]
        xyz_new[i,1] = xyz[i,1]+wrapy[i]
        xyz_new[i,2] = xyz[i,2]+wrapz[i]
        xmin=-Lx+(xz-xy*yz)*xyz_new[i,2]+xy*xyz_new[i,1]
        xmax= Lx+(xz-xy*yz)*xyz_new[i,2]+xy*xyz_new[i,1]
        ymin=-Ly+yz*xyz_new[i,2]
        ymax= Ly+yz*xyz_new[i,2]
        while ((xyz[i,0]+wrapx[i]) < xmin):
#            print('1',xyz[i])
            wrapx[i] -= a[0]
            xyz_wp[i,0] += 1
            if ((xyz[i,0]+wrapx[i]) >= xmin):
                break
        while ( (xyz[i,0]+wrapx[i]) > xmax):
#            print('2',xyz[i])
            wrapx[i] += a[0]
            xyz_wp[i,0] -= 1
            if ((xyz[i,0]+wrapx[i]) <= xmax):
                break
        while ((xyz[i,1]+wrapy[i])  < ymin):
#            print('3',xyz[i])
            wrapy[i] += b[1]
            xyz_wp[i,1] -= 1
            if ((xyz[i,1]+wrapy[i]) >= ymin):
                break
        while ((xyz[i,1]+wrapy[i]) > ymax):
#            print('4',xyz[i])
            wrapy[i] -= b[1]
            xyz_wp[i,1] += 1
            if ( (xyz[i,1]+wrapy[i]) <= ymax):
                break
        while ((xyz[i,2]+wrapz[i]) < zmin):
#            print('5',xyz[i])
            wrapx[i] += c[0]
            wrapz[i] += c[2]
            xyz_wp[i,2] -= 1
            if ((xyz[i,2]+wrapz[i]) >= zmin):
                break
        while ((xyz[i,2]+wrapz[i]) > zmax):
#            print('6',xyz[i])
            wrapx[i] -= c[0]
            wrapz[i] -= c[2]
            xyz_wp[i,0] += 1
            if ((xyz[i,2]+wrapz[i]) <= zmax):
                break
for i in range(l2_num):
    aax[i] = xyz_new[i,0] + T2[0]           
    aay[i] = xyz_new[i,1] + T2[1]    
    aaz[i] = xyz_new[i,2] + T2[2]    
write_xyz('mirror.xyz',atype,aax,aay,aaz)
ax=[];ay=[];az=[];at=[]
for i in range(l2_num):
    ax.append(xyz_new[i,0] + T2[0])           
    ay.append(xyz_new[i,1] + T2[1])    
    az.append(xyz_new[i,2] + T2[2])    
    at.append(atype[i])
for i in range(l2_num):
    ax.append(xyz_new[i,0] + T2[0]+a[0])           
    ay.append(xyz_new[i,1] + T2[1]+a[1])    
    az.append(xyz_new[i,2] + T2[2]+a[2])    
    at.append(atype[i])
for i in range(l2_num):
    ax.append(xyz_new[i,0] + T2[0]+b[0])           
    ay.append(xyz_new[i,1] + T2[1]+b[1])    
    az.append(xyz_new[i,2] + T2[2]+b[2])    
    at.append(atype[i])
for i in range(l2_num):
    ax.append(xyz_new[i,0] + T2[0]+a[0]+b[0])           
    ay.append(xyz_new[i,1] + T2[1]+a[1]+b[1])    
    az.append(xyz_new[i,2] + T2[2]+a[2]+b[2])    
    at.append(atype[i])
write_xyz('mirror2.xyz',at,ax,ay,az)
print('fini')