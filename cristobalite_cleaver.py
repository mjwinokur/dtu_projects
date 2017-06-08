#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 11:42:46 2017

@author: winokur
"""#
#from mw_library import q_mult
#from mw_library import normalize
#from mw_library import q_conjugate
#from mw_library import axisangle_to_q
#from mw_library import q_to_axisangle
#from mw_library import qv_mult
from mw_library import read_xyz
from mw_library import read_txyz
from mw_library import write_xyz
from mw_library import write_txyz
from mw_library import write_tinker_key
#from mw_library import write_txyz
#from mw_library import iofiles
#from rbuildTx import rebuild_Tx
import numpy as np
#import sys
#import random
import os
#import scipy.optimize
#import functools
#import math
#import os
#from mayavi.mlab import *
#########################################################################
#
# Main:  This program builds a cristobalite supercell from a base unit cell
# in Tinker format 
#
#########################################################################
#xyzname = "cristobalite.xyz" # Base
input_name='ctest3.xyz_6'
abt=[];astring=[]
l2_num,aax,aay,aaz,atype,astring,abt,uc,header = read_txyz(input_name)
#It seems the best way is to identify pairs which are along a-axis
#Identify Si atoms
aaxis = np.array([1.0,0.0,0.0])
baxis = np.array([0.0,1.0,0.0])
caxis = np.array([0.0,0.0,1.0])
atest = uc[0]-2.0  # minumum distance for breaking bond
si_site = []
break_n = []
br_atm = []
for i in range(l2_num):
    if (abt[i]==19):
        si_site.append(i)
# Now calculate pair distances
for i in si_site:
    xyz0 = np.array([aax[i],aay[i],aaz[i]])
    newline = ' '.join(astring[i].split())
    mylist = newline.split(" ")
    nrem = 0
    mylistp = []
    for j in mylist:
        mylistp.append(int(j))
    for j in mylistp:
        k = j-1
        xyz1 = np.array([aax[k],aay[k],aaz[k]])
        delta=xyz1-xyz0
        dist = np.sqrt(np.dot(delta,delta))
        if (dist > 2.0):
            proj = np.dot(aaxis,delta)
#            print 'b',i,k,atest,proj
#            raw_input()
            if (proj > atest):
                break_n.append(k)
                br_atm.append(i)
                mylist.remove(str(k+1))
                nrem += 1
    if (nrem > 0):  # Need to revise astring element
        newstring = ""
        for k in mylist:
            tstring = "%5s" % (k) 
            newstring = newstring + tstring
        astring[i] = newstring
#            print i,k,dist,proj
#        print 'nrem',nrem
#        print i,mylist,astring[i]
# Now to remove the pair from the oyxgen atom list
oxy_atm = []
for j in range(len(break_n)):
    i = break_n[j]
    newline = ' '.join(astring[i].split())
    mylist = newline.split(" ")
    mylist.remove(str(br_atm[j]+1))
    newstring = ""
    for k in mylist:
        tstring = "%5s" % (k) 
        newstring = newstring + tstring
    astring[i] = newstring
# Hack to deal with np array limitations
ax=[];ay=[];az=[];
for i in range(l2_num):
    ax.append(aax[i])
    ay.append(aay[i])
    az.append(aaz[i])
# Now to protonate the oxygens
for i in break_n:
    if (abt[i] == 6): # oxygen
        l2_num += 1
        atype.append('H')
        abt.append(21)  # Alcohol
        ax.append(aax[i]+1.)
        ay.append(aay[i])
        az.append(aaz[i])
        tstring = "%4s" % (i+1) 
        astring.append(tstring)
        tstring = "%5s" % (l2_num) 
        astring[i] += tstring
# Now to protonate the silicons
hs=0.7
for i in br_atm:
    if (abt[i] == 19): # silicons
        l2_num += 1
        atype.append('H')
        abt.append(5)  # Silane
        ax.append(aax[i]-0.7)
        ay.append(aay[i]+hs)
        hs *= -1.
        az.append(aaz[i])
        tstring = "%4s" % (i+1) 
        astring.append(tstring)
        tstring = "%5s" % (l2_num) 
        astring[i] += tstring
uc[0]=uc[0]+20.    
write_txyz('ctest4.txyz',l2_num,atype,astring,ax,ay,az,abt,uc,header)
write_tinker_key('ctest4.key',uc)
print 'Done'