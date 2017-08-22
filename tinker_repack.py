#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 13:38:56 2017

@author: winokur
"""
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 13:00:37 2017

@author: winokur
"""
#
#from mw_library import q_mult
#from mw_library import normalize
#from mw_library import q_conjugate
#from mw_library import axisangle_to_q
#from mw_library import q_to_axisangle
#from mw_library import qv_mult
from mw_library import read_txyz
from mw_library import write_txyz
from mw_library import iofiles
#from mw_library import merge_txyz
#from rbuildTx import rebuild_Tx
import numpy as np
import sys
#import random
#import scipy.optimize
#import functools
#import os
#from mayavi.mlab import *
#########################################################################
#
# Main:  Map NaT2 or NaT3 molecules contiguously and shrink the unit cell along the a-axis

#########################################################################
#
# open base momomer as an input 
txyzname = "S1A_101_15n/mol_testa.xyz_3" # Base monomer in tinker format and includes the unit cell properties on line 2
#
# Two output files,  base name (to the left must be the same)
#
# 
txyznewname = "S1A_101_15n/mol_testa.xyz_4"  # or just ztest.xyz   ########## Look at me!
tnewkey="S1A_101_15n/mol_testa.key"  ######### Look at me!
# Repeat of the unit cell along the principle lattice vectors,  currently 3x3x3
#alpha=0.; beta=0.; gamma = 0.0
#iflip = 'all'
monomer = 'T2'
if (len(sys.argv[1:])== 0):
    print "Using defaults", txyzname,txyznewname,tnewkey
#    print "With a cell that is ",lenx,"x",leny,"x",lenz
    print "And with no flips or custom build options"
else: 
    ttxyzname,ttxyznewname,text = iofiles(sys.argv[1:])
    if (ttxyznewname != ""):
        txyznewname = ttxyznewname+'.txyz'
        tnewkey = ttxyznewname+'.key'
    if (ttxyzname !=""):
        txyzname = ttxyzname
    print "Files are", txyzname,txyznewname,tnewkey
    if (text != ""):
        print "with special:",text
        alist = text.split()
        j = 0
 #       print 'alist:',alist,j
        for i in alist:
            if (alist[j]== 'option'):
                iflip ='option'
            j += 1
#
l2_num,aax,aay,aaz,atype,astring,abt,uc,header = read_txyz(txyzname)
cosa=np.cos(np.radians(uc[3])) #; sina=np.sin(np.radians(uc[3]))
cosb=np.cos(np.radians(uc[4])) #; sinb=np.sin(np.radians(uc[4]))
cosg=np.cos(np.radians(uc[5])); sing=np.sin(np.radians(uc[5]))
Vdivab=uc[2]*np.sqrt(1.-cosa*cosa-cosb*cosb-cosg*cosg+2.*cosa*cosb*cosg)
a =[uc[0],0.,0.]
b =[uc[1]*cosg,uc[1]*sing,0.]
c =[(uc[2]*cosb),uc[2]*(cosa-cosb*cosg)/sing,Vdivab/(sing)]
# 
# find the center of each molecule
#
natm=96
if (monomer == 'T3'):
    natm=110
m2ct=l2_num/natm
for i in range(m2ct):
    datamean = []
    i0=i*natm
    if (monomer == 'T2'):
        xt = []; yt = []; zt = []
        i1=i0+24
        i2=i0+48
        i3=i0+72
        i4=i0+96
        xt = np.append(np.asarray(aax[i0:i1],dtype=float),np.asarray(aax[i2:i3],dtype=float))
        yt = np.append(np.asarray(aay[i0:i1],dtype=float),np.asarray(aay[i2:i3],dtype=float))
        zt = np.append(np.asarray(aaz[i0:i1],dtype=float),np.asarray(aaz[i2:i3],dtype=float))
        data = np.concatenate((xt[:,np.newaxis],yt[:,np.newaxis],zt[:,np.newaxis]),axis=1)
# Calculate the mean of the points, i.e. the 'center' of the cloud
        [xm,ym,zm]=data.mean(axis=0)
        if (xm <0.0):
            aax[i0:i1] += a[0] 
            aay[i0:i1] += a[1] 
            aaz[i0:i1] += a[2] 
            aax[i2:i3] += a[0] 
            aay[i2:i3] += a[1] 
            aaz[i2:i3] += a[2] 
#
        xt = []; yt = []; zt = []
        xt = np.append(np.asarray(aax[i1:i2],dtype=float),np.asarray(aax[i3:i4],dtype=float))
        yt = np.append(np.asarray(aay[i1:i2],dtype=float),np.asarray(aay[i3:i4],dtype=float))
        zt = np.append(np.asarray(aaz[i1:i2],dtype=float),np.asarray(aaz[i3:i4],dtype=float))
        data = np.concatenate((xt[:,np.newaxis],yt[:,np.newaxis],zt[:,np.newaxis]),axis=1)
# Calculate the mean of the points, i.e. the 'center' of the cloud
        [xm,ym,zm]=data.mean(axis=0)
        if (xm <0.0):
            aax[i1:i2] += a[0] 
            aay[i1:i2] += a[1] 
            aaz[i1:i2] += a[2] 
            aax[i3:i4] += a[0] 
            aay[i3:i4] += a[1] 
            aaz[i3:i4] += a[2] 
    elif (monomer == 'T3'):
        i1=i0+55
        i2=i0+110
        xt = np.asarray(aax[i0:i1],dtype=float)
        yt = np.asarray(aay[i0:i1],dtype=float)
        zt = np.asarray(aaz[i0:i1],dtype=float)
        data = np.concatenate((xt[:,np.newaxis],yt[:,np.newaxis],zt[:,np.newaxis]),axis=1)
# Calculate the mean of the points, i.e. the 'center' of the cloud
        [xm,ym,zm]=data.mean(axis=0)
        if (xm <0.0):
            aax[i0:i1] += a[0] 
            aay[i0:i1] += a[1] 
            aaz[i0:i1] += a[2] 
#
        xt = np.asarray(aax[i1:i2],dtype=float)
        yt = np.asarray(aay[i1:i2],dtype=float)
        zt = np.asarray(aaz[i1:i2],dtype=float)
        data = np.concatenate((xt[:,np.newaxis],yt[:,np.newaxis],zt[:,np.newaxis]),axis=1)
# Calculate the mean of the points, i.e. the 'center' of the cloud
        [xm,ym,zm]=data.mean(axis=0)
        if (xm <0.0):
            aax[i0:i1] += a[0] 
            aay[i0:i1] += a[1] 
            aaz[i0:i1] += a[2] 

#    print i, temp,abt[i]
#    raw_input()
##########################
# Now to save the positions   
##########################
#
write_txyz(txyznewname,l2_num,atype,astring,aax,aay,aaz,abt,uc,header)
#
#
f = open(tnewkey,'w')
#f.write('# Force Field Selection \nPARAMETERS        /home/winokur/MolecularTools/ffe/../tinker/params/mm3.prm \n')
f.write('# Force Field Selection \nPARAMETERS        /home/winokur/MolecularTools/ffe/../tinker/params/oplsaa.prm \n')
f.write('# Crystal Lattice And Periodic Boundary \n\n')
f.write('A-AXIS   '+str(uc[0])+' \n')
f.write('B-AXIS   '+str(uc[1])+' \n')
f.write('C-AXIS   '+str(uc[2])+' \n')
f.write('ALPHA    '+str(uc[3])+' \n')
f.write('BETA     '+str(uc[4])+' \n')
f.write('GAMMA    '+str(uc[5])+' \n')
f.write('\n')
#for i in range(nspecial):
#    f.write('RESTRAIN-POSITION '+ str(nsp[i]))
#    f.write('#PISYSTEM '+ mynewstring+'\n')
f.close()  # close write file
# fini

