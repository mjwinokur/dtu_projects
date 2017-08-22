#
#from mw_library import q_mult
#from mw_library import normalize
#from mw_library import q_conjugate
from mw_library import axisangle_to_q
#from mw_library import q_to_axisangle
from mw_library import qv_mult
from mw_library import read_txyz
from mw_library import write_txyz
import numpy as np
#import scipy.optimize
#import functools
#import math
#import os
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
#
txyzname='NaT2mw.txyz'
l2_num,aax,aay,aaz,atype,astring,abt,uc,header = read_txyz(txyzname)
#aax = np.asarray(aax,dtype=float)
#aay = np.asarray(aay,dtype=float)
#aaz = np.asarray(aaz,dtype=float)
tilta = 90.; tiltb=0.; tiltc=0.  # 90 degrees about the a-axis
ax,ay,az = reorient_one(tilta,tiltb,tiltc,uc,aax,aay,aaz)
write_txyz('testing.txyz_1',l2_num,atype,astring,ax,ay,az,abt,uc,header)
tilta =  0.; tiltb=90.; tiltc=0. # 90 degrees about the b-axis
ax,ay,az = reorient_one(tilta,tiltb,tiltc,uc,aax,aay,aaz)
write_txyz('testing.txyz_2',l2_num,atype,astring,ax,ay,az,abt,uc,header)
tilta =  0.; tiltb=0.; tiltc=90. # 90 degreees about the c-axis
ax,ay,az = reorient_one(tilta,tiltb,tiltc,uc,aax,aay,aaz)
write_txyz('testing.txyz_3',l2_num,atype,astring,ax,ay,az,abt,uc,header)
print 'Done'
#