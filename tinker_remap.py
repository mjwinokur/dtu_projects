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
#import math
#import os
#from mayavi.mlab import *
#########################################################################
#
# Main:  Convert mm3 atom referencing to oplsaa
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
# open base momomer as an input 
#txyzname = "NaT3_s_tinker_min.xyz" # From tinker minimization
#txyzname = "NaT3_s_tinker_min.xyz_2" # From tinker minimization
txyzname = "NaT2101a.xyz_17" # Base monomer in tinker format and includes the unit cell properties on line 2
#
# Two output files,  base name (to the left must be the same)
#
# 
txyznewname = "mol_testa.xyz"  # or just ztest.xyz   ########## Look at me!
tnewkey="mol_testa.key"  ######### Look at me!
# Repeat of the unit cell along the principle lattice vectors,  currently 3x3x3
#alpha=0.; beta=0.; gamma = 0.0
#iflip = 'all'
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
# 
atomout =np.zeros(1000,dtype=int)
"""
atomout[2]= 2 # 90   48 Aromatic C
atomout[5]= 5 # 91   49 Aromatic H
atomout[42]= 42 #26   16  Thioether
atomout[6]=  6 # 870   
atomout[21]= 21 # 870   45 Alkyl Silane H-C-Si
atomout[19]= 19 #868  108 Alkyl Silane R2SiH2
"""
atomout[2]  =  90 # 90   48 Aromatic C
atomout[5]  =  91 # 91   49 Aromatic H
atomout[42] =  26 #26   16  Thioether
atomout[6]  = 908 # Oxygen  
atomout[19] = 907 # Si
#atomout[6]=  520 # Oxygen  
atomout[21] =  85 # 870   45 Alkyl Silane H-C-Si
#atomout[19]= 868 #868  108 Alkyl Silane R2SiH2
for i in range(l2_num):
    temp = abt[i]
    abt[i]=atomout[temp]
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
