#
#from mw_library import q_mult
#from mw_library import normalize
#from mw_library import q_conjugate
from mw_library import axisangle_to_q
#from mw_library import q_to_axisangle
from mw_library import qv_mult
#from mw_library import read_txyz
#from mw_library import write_txyz
import numpy as np
#import scipy.optimize
#import functools
#import math
#import os
#from mayavi.mlab import *
#########################################################################
#
# Main:  This program builds a NaT2 supercell from a base unit cell
# in Tinker format and then, using hoomd_gsd_test, takes this information and 
# that from openbabel to code two files for use in hoomd.  The first files, a 
# a gsd format file
# contains the information for the atoms and connectivity and encodes for
# periodic boundary conditions.  The second, a Python dump, contains the 
# interaction coefficients and other key details.   
#
#########################################################################
def plane(x, y, params):
    a = params[0]
    b = params[1]
    c = params[2]
    z = a*x + b*y + c
    return z

def error(params, points):
    result = 0
    for (x,y,z) in points:
        plane_z = plane(x, y, params)
        diff = abs(plane_z - z)
        result += diff**2
    return result
#
def rebuild_Tx(tilta,tiltb,tiltc,uc,aax,aay,aaz,l2_num):
# uc is the conventional way that a 3D lattice is stated: a,b,c,alpha,beta,gamma in crystallography
# a is always paralel to the x-axis
# b us always in the xy plane
# tilta, tiltb, tilc are the amount of tilt (in degrees) around the a, b and c axes respectively
# aax,aay,aaz are 1D arrays of the x,y and z coordinates of the atoms
# l2_num is the number of atoms in a molecule
# As it is now configured there are explicit assumptions for the NaT2 and NaT3 monomer construction
# In particular there are two monomers per unit cell
# Currently the approach is to find the xyz geometric center of each molecular unit with the resulting value appended to datamean
# So, for the first molecular of the pair, datamean[0] = [x-center,y-center,z-center]
# and for the second molecular of the pair, datamean[1] = [x-center,y-center,z-center]
# It also calculated the best-fit line that passes through all the atoms (i.e., the molecular axis) but this 
# result bfline is not used further on as that aspect is currently commented out
# 
#
#nct=lenx*leny*lenz*l2_num
    n_monomers = len(aax)/l2_num
    cosa=np.cos(np.radians(uc[3])) #; sina=np.sin(np.radians(uc[3]))
    cosb=np.cos(np.radians(uc[4])) #; sinb=np.sin(np.radians(uc[4]))
    cosg=np.cos(np.radians(uc[5])); sing=np.sin(np.radians(uc[5]))
    Vdivab=uc[2]*np.sqrt(1.-cosa*cosa-cosb*cosb-cosg*cosg+2.*cosa*cosb*cosg)
    a =[uc[0],0.,0.]
    b =[uc[1]*cosg,uc[1]*sing,0.]
    c =[(uc[2]*cosb),uc[2]*(cosa-cosb*cosg)/sing,Vdivab/(sing)]
#c =[(uc[2]*np.cos(np.radians(uc[4]))),0.,(uc[2]*np.sin(np.radians(uc[4])))]
    q1R1 = axisangle_to_q(a,np.deg2rad(tilta))
    q1R2 = axisangle_to_q(b,np.deg2rad(tiltb))
    q1R3 = axisangle_to_q(c,np.deg2rad(tiltc))
    itt = l2_num/2
    for i in range (n_monomers):
        datamean = []
        bfline = []
        if (l2_num == 96):
            i1 = i*l2_num
            i2 = i1 +24
            i3 = i1 +48
            i4 = i1 +72
            i5 = i1 +96
            xt = np.append(np.asarray(aax[i1:i2],dtype=float),np.asarray(aax[i3:i4],dtype=float))
            yt = np.append(np.asarray(aay[i1:i2],dtype=float),np.asarray(aay[i3:i4],dtype=float))
            zt = np.append(np.asarray(aaz[i1:i2],dtype=float),np.asarray(aaz[i3:i4],dtype=float))
            data = np.concatenate((xt[:,np.newaxis],yt[:,np.newaxis],zt[:,np.newaxis]),axis=1)
        # Calculate the mean of the points, i.e. the 'center' of the cloud
            datamean.append(data.mean(axis=0))
# Do an SVD on the mean-centered data.
            uu, dd, vv = np.linalg.svd(data - datamean[0])
            bfline.append(vv[0])
#
            xt = np.append(np.asarray(aax[i2:i3],dtype=float),np.asarray(aax[i4:i5],dtype=float))
            yt = np.append(np.asarray(aay[i2:i3],dtype=float),np.asarray(aay[i4:i5],dtype=float))
            zt = np.append(np.asarray(aaz[i2:i3],dtype=float),np.asarray(aaz[i4:i5],dtype=float))
            data = np.concatenate((xt[:,np.newaxis],yt[:,np.newaxis],zt[:,np.newaxis]),axis=1)
# Calculate the mean of the points, i.e. the 'center' of the cloud
            datamean.append(data.mean(axis=0))
# Do an SVD on the mean-centered data.
            data = data - datamean[0]
            uu, dd, vv = np.linalg.svd(data)
            bfline.append(vv[0])
# vv[0] is the best-fit line
            imap = np.empty(l2_num,dtype=int)
            for i in range(24):  # the atoms of a monomer are not contiguous
                imap[i]=i+i1
                imap[i+24]=i+itt+i1
                imap[i+48]=i+24+i1
                imap[i+72]=i+itt+24+i1
        if (l2_num == 110):
            i1 = i*l2_num
            i2 = i1 +55
            i3 = i1 +110
            xt = np.asarray(aax[i1:i2],dtype=float)
            yt = np.asarray(aay[i1:i2],dtype=float)
            zt = np.asarray(aaz[i1:i2],dtype=float)
            data = np.concatenate((xt[:,np.newaxis],yt[:,np.newaxis],zt[:,np.newaxis]),axis=1)
# Calculate the mean of the points, i.e. the 'center' of the cloud
            datamean.append(data.mean(axis=0))
# Do an SVD on the mean-centered data.
            uu, dd, vv = np.linalg.svd(data - datamean[0])
            bfline.append(vv[0])
#
            xt = np.asarray(aax[i2:i3],dtype=float)
            yt = np.asarray(aay[i2:i3],dtype=float)
            zt = np.asarray(aaz[i2:i3],dtype=float)
            data = np.concatenate((xt[:,np.newaxis],yt[:,np.newaxis],zt[:,np.newaxis]),axis=1)
# Calculate the mean of the points, i.e. the 'center' of the cloud
            datamean.append(data.mean(axis=0))
# Do an SVD on the mean-centered data.
            data = data - datamean[0]
            uu, dd, vv = np.linalg.svd(data)
            bfline.append(vv[0])  # This is the best-fit line
            imap = np.empty(l2_num,dtype=int)
            for i in range(l2_num):  # the atoms of a monomer are contiguous
                imap[i]=i+i1
#
# Vectors perpendicular to the bfline
# Find best-fit plane
#fun = functools.partial(error, points=data)
#params0 = [0, 0, 0]
#res = scipy.optimize.minimize(fun, params0)
#a = res.x[0]; b = res.x[1]; c = res.x[2]
#point  = np.array([0.0, 0.0, c])
#norm_out = np.cross([1,0,a], [0,1,b])  #Normal to the plane
#d = -point.dot(norm_out)
###########################  Commented out   
# The oligomer orientation may need to be altered which uses bfline
#
# Use quaternions to rotate molecules
#  With zero this does nothing but for NaT3 it will matter
#q1 = axisangle_to_q([aax[i1]-aax[i2],aay[i1]-aay[i2],aaz[i1]-aaz[i2]],0.0)
#        q1 = axisangle_to_q(bfline[0],0.0)
# This will rotate one of the two monomer 180 deg to give head to head
#q1R = axisangle_to_q([aax[i1]-aax[i2],aay[i1]-aay[i2],aaz[i1]-aaz[i2]],np.pi)
#        q1R = axisangle_to_q(bfline[0],np.pi)
#########################  Comment out
        [cx1,cy1,cz1] = datamean[0]
        [cx2,cy2,cz2] = datamean[1]
        for j in range(itt):
            i = imap[j] 
            temp = [aax[i]-cx1,aay[i]-cy1,aaz[i]-cz1]
            temp = qv_mult(q1R1,(temp[0],temp[1],temp[2])) 
            temp = qv_mult(q1R2,(temp[0],temp[1],temp[2])) 
            temp = qv_mult(q1R3,(temp[0],temp[1],temp[2])) 
            aax[i]=temp[0]+cx1;aay[i]=temp[1]+cy1;aaz[i]=temp[2]+cz1;        
#
            i = imap[j+itt] 
            temp = [aax[i]-cx2,aay[i]-cy2,aaz[i]-cz2]
            temp = qv_mult(q1R1,(temp[0],temp[1],temp[2])) 
            temp = qv_mult(q1R2,(temp[0],temp[1],temp[2])) 
            temp = qv_mult(q1R3,(temp[0],temp[1],temp[2])) 
            aax[i]=temp[0]+cx2;aay[i]=temp[1]+cy2;aaz[i]=temp[2]+cz2;
#   
#    print i, xx[i],yy[i],zz[i]
    return aax,aay,aaz
#
def reorient_one(tilta,tiltb,tiltc,uc,aax,aay,aaz):
#
# uc is the conventional way that a 3D lattice is stated: a,b,c,alpha,beta,gamma in crystallography
# a is always paralel to the x-axis
# b us always in the xy plane
# tilta, tiltb, tilc are the amount of tilt (in degrees) around the a, b and c axes respectively
# aax,aay,aaz are 1D arrays of the x,y and z coordinates of the atoms
# Here the approach is to find the xyz geometric center of entire input array
# 
#
    cosa=np.cos(np.radians(uc[3])) #; sina=np.sin(np.radians(uc[3]))
    cosb=np.cos(np.radians(uc[4])) #; sinb=np.sin(np.radians(uc[4]))
    cosg=np.cos(np.radians(uc[5])); sing=np.sin(np.radians(uc[5]))
    Vdivab=uc[2]*np.sqrt(1.-cosa*cosa-cosb*cosb-cosg*cosg+2.*cosa*cosb*cosg)
    a =[uc[0],0.,0.]
    b =[uc[1]*cosg,uc[1]*sing,0.]
    c =[(uc[2]*cosb),uc[2]*(cosa-cosb*cosg)/sing,Vdivab/(sing)]
#c =[(uc[2]*np.cos(np.radians(uc[4]))),0.,(uc[2]*np.sin(np.radians(uc[4])))]
    q1R1 = axisangle_to_q(a,np.deg2rad(tilta))
    q1R2 = axisangle_to_q(b,np.deg2rad(tiltb))
    q1R3 = axisangle_to_q(c,np.deg2rad(tiltc))
    cx1=np.mean(aax)
    cy1=np.mean(aay)
    cz1=np.mean(aaz)
    for i in range(len(aax)):
            temp = [aax[i]-cx1,aay[i]-cy1,aaz[i]-cz1]
            temp = qv_mult(q1R1,(temp[0],temp[1],temp[2])) 
            temp = qv_mult(q1R2,(temp[0],temp[1],temp[2])) 
            temp = qv_mult(q1R3,(temp[0],temp[1],temp[2])) 
            aax[i]=temp[0]+cx1;aay[i]=temp[1]+cy1;aaz[i]=temp[2]+cz1;        
#    print i, xx[i],yy[i],zz[i]
    return aax,aay,aaz
#
# Test code
"""
txyzname='NaT2mw.txyz'
l2_num,aax,aay,aaz,atype,astring,abt,uc,header = read_txyz(txyzname)
#aax = np.asarray(aax,dtype=float)
#aay = np.asarray(aay,dtype=float)
#aaz = np.asarray(aaz,dtype=float)
tilta = 90.; tiltb=0.; tiltc=0.
ax,ay,az = reorient_one(tilta,tiltb,tiltc,uc,aax,aay,aaz)
write_txyz('testing.txyz_1',l2_num,atype,astring,aax,aay,aaz,abt,uc,header)
tilta =  0.; tiltb=90.; tiltc=0.
ax,ay,az = reorient_one(tilta,tiltb,tiltc,uc,aax,aay,aaz)
write_txyz('testing.txyz_2',l2_num,atype,astring,aax,aay,aaz,abt,uc,header)
tilta =  0.; tiltb=0.; tiltc=90.
ax,ay,az = reorient_one(tilta,tiltb,tiltc,uc,aax,aay,aaz)
write_txyz('testing.txyz_3',l2_num,atype,astring,aax,aay,aaz,abt,uc,header)
print 'Done'
"""
#