#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 11:54:24 2017

@author: winokur
"""

#
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
xyzname = "9008234.xyz" # Base
# Repeat of the unit cell along the principle lattice vectors,  currently 3x3x3
lenx=2; leny=3; lenz=3
a = 7.1666; b = 7.1666; c = 7.1666
#a = 5.067; b = 5.067; c = 5.067
alpha=90.; beta=90.; gamma = 90.0
#
l_num,ax,ay,az,at = read_xyz(xyzname)
a2x=[];a2y=[];a2z=[];a2t=[]
# The cif file has multiple fractional atoms in a site.  We just need to keep an average
use = np.zeros(l_num)
for i in range(l_num):
    if (use[i] == 0):
        xl=[];yl=[];zl=[]
#        use[i]=1
        xl.append(ax[i]);yl.append(ay[i]);zl.append(az[i]);
        x0=xl[0];y0=yl[0];z0=zl[0]    
        for j in range(i+1,l_num):
            if (use[j]== 0):
                x1 = x0-ax[j]
                y1 = y0-ay[j]
                z1 = z0-az[j]
                dist = np.sqrt(x1*x1+y1*y1+z1*z1)
                if (dist < 0.9):
                    xl.append(ax[j])
                    yl.append(ay[j])
                    zl.append(az[j])
                    use[j]=1
        a2x.append(np.mean(xl)) 
        a2y.append(np.mean(yl)) 
        a2z.append(np.mean(zl)) 
        a2t.append(at[i])
#
a2x =np.asarray(a2x)
a2y =np.asarray(a2y)
a2z =np.asarray(a2z)
x0 = np.mean(a2x);x1=x0-a*0.5 
y0 = np.mean(a2y);y1=y0-b*0.5 
z0 = np.mean(a2z);z1=z0-c*0.5 
#recenter atoms and keep those in unit cell
l_num=len(a2x)
aax=[];aay=[];aaz=[];atype=[]
for i in range(l_num):
#    a2x[i] -= x1; a2y[i] -= y1;a2z[i] -= z1
    x0 = a2x[i]-x1;y0 = a2y[i]-y1; z0 =a2z[i]-z1
    if ((x0 >= 0.00 and x0 < a) and (y0 >= 0.00 and y0 < b) and (z0 >= 0.00 and z0 < c)):
        aax.append(x0)
        aay.append(y0)
        aaz.append(z0)
        atype.append(a2t[i])
#write_xyz('ctest.xyz',atype,aax,aay,aaz)
#
# eliminate duplicates 
#
l_num=len(aax)
keep = np.zeros(l_num)
for i in range(l_num):    
    x0=aax[i];y0=aay[i];z0=aaz[i]
    if (x0 > a*0.5):
        x0 -= a
    if (y0 > b*0.5):
        y0 -= b
    if (z0 > c*0.5):
        z0 -= c
    for j in range(i+1,l_num):
        x1=aax[j];y1=aay[j];z1=aaz[j]
        if (x1 > a*0.5):
            x1 -= a
        if (y1 > b*0.5):
            y1 -= b
        if (z1 > c*0.5):
            z1 -= c
        x2 = x1-x0
        y2 = y1-y0
        z2 = z1-z0
        dist = np.sqrt(x2*x2+y2*y2+z2*z2)
        if (dist < 0.05):
#            print i,j,'duplicate'
            keep[j]=1
a2x=[];a2y=[];a2z=[];a2t=[]
for i in range(l_num):
    if (keep[i] ==0):
        a2x.append(aax[i])
        a2y.append(aay[i])
        a2z.append(aaz[i])
        a2t.append(atype[i])
l_num=len(a2x)
#
# Replicate the unit cell
#
aax=[];aay=[];aaz=[];atype=[]
for i in range(lenx+2):
    x0=a*float(i-1)
    for j in range(leny+2):
        y0=b*float(j-1)
        for k in range (lenz+2):
            z0=c*float(k-1)
            for ii in range(l_num):
                x=x0+a2x[ii]
                y=y0+a2y[ii]
                z=z0+a2z[ii]
                aax.append(x)
                aay.append(y)
                aaz.append(z)
                atype.append(a2t[ii])
#
# Reorder sequence so that those atom outside central range are sequenced
# 1: positions less that zero
# 2: keepers in replicated cell
# 3: positions beyond replicated cell 
l2_num=len(aax)
a0seq=[];a1seq=[];a2seq=[]
amax=a*float(lenx);bmax=b*float(leny);cmax=c*float(lenz)
for i in range(l2_num):
    if (aax[i] < 0.0 or aay[i] < 0.0 or aaz[i] < 0.0 ):
        a0seq.append(i)
    elif (aax[i] >= amax or aay[i] >= bmax or aaz[i] >= cmax ):
        a2seq.append(i)
    else:
        a1seq.append(i)
a2x=[];a2y=[];a2z=[];a2t=[]
for i in a0seq:
    a2x.append(aax[i]);a2y.append(aay[i]);a2z.append(aaz[i]);a2t.append(atype[i])
for i in a1seq:
    a2x.append(aax[i]);a2y.append(aay[i]);a2z.append(aaz[i]);a2t.append(atype[i])
for i in a2seq:
    a2x.append(aax[i]);a2y.append(aay[i]);a2z.append(aaz[i]);a2t.append(atype[i])
write_xyz('ctest.xyz',a2t,a2x,a2y,a2z)
#write_xyz('ctest.xyz',atype,aax,aay,aaz)
command='obabel -ixyz ctest.xyz -otxyz -Octest.txyz'
print ' cmd:',command
os.system(command)
# obabel is not perfect and breaks the txyz file
# We need to fix it and put in proper periodic boundary conditions  
l2_num,aax,aay,aaz,atype,astring,abt,uc,header = read_txyz('ctest.txyz')
# now to fix the mistakes
nmiss=0
for i in range(l2_num):
    if (abt[i] == 0 ):        
        abt[i] = 19 # Silicons are set to zero instead of mm3
    if (astring[i]== -1):
        astring[i]=' '
        nmiss += 1
uc = np.zeros(6)
uc[0]=a*float(lenx+2)+20.;uc[1]=b*float(leny+2)+20.;uc[2]=c*float(lenz+2)+20.;
uc[3]=90.0;uc[4]=90.0;uc[5]=90.0
#
# Add missing nearest neighbors...if any; current version seems to avoid this
if (nmiss > 0):
    for i in range(l2_num):        
        if (astring[i]== ' '):
            for j in range(l2_num):
                if (abt[j]==19 and astring[i]==' '):
                    x0 = aax[j]-aax[i]
                    y0 = aay[j]-aay[i]
                    z0 = aaz[j]-aaz[i]
                    dist=np.sqrt(x0*x0+y0*y0+z0*z0)
                    if (dist < 1.7):
#                        print i,j,dist,atype[j]
                        astring[i]="%4s" % (str(j+1))
l2_num=len(aax)
write_txyz('ctest2.txyz',l2_num,atype,astring,aax,aay,aaz,abt,uc,header)
write_tinker_key('ctest2.key',uc)
# Now we need to reindex nearest neighbors consistent with periodic boundary conditions
a2x=[];a2y=[];a2z=[];a2t=[];ab2t=[];a2string=[]
pair1=[];pair2=[]
uc[0]=a*float(lenx);uc[1]=b*float(leny);uc[2]=c*float(lenz);
j = len(a0seq)
jj = len(a1seq)
k = j+ jj
for i in range(j,k):
    a2x.append(aax[i]);a2y.append(aay[i]);a2z.append(aaz[i]);ab2t.append(abt[i])
# decode 
testct = 0
for i in range(j,k):
    newline = ' '.join(astring[i].split())
    mylist = newline.split(" ")
    alen = len(mylist)
    newstring =''        
    for n in range(alen):
        aval = int(mylist[n])-1 # Index shift for array
        bval = aval-j
        if ((bval < 0 or bval > jj) and abt[i] == 19 ): # need to find proper nearest neighbor
#            print i,i-j+1, mylist
            x0 = aax[aval]
            y0 = aay[aval]
            z0 = aaz[aval]
            if (x0 < 0.0):
                x0 += uc[0]
            if (y0 < 0.0):
                y0 += uc[1]
            if (z0 < 0.0):
                z0 += uc[2]
            if (x0 >= uc[0]):
                x0 -= uc[0]
            if (y0 >= uc[1]):
                y0 -= uc[1]
            if (z0 >= uc[2]):
                z0 -= uc[2]
            for kk in range(jj):
                x1 = a2x[kk] - x0
                y1 = a2y[kk] - y0
                z1 = a2z[kk] - z0
                dist = np.sqrt(x1*x1+y1*y1+z1*z1)
                if (dist < 0.01): 
                    pair1.append(i-j)
                    pair2.append(kk)
#                    print i-j+1,kk+1
        elif (bval >=0 and bval <= jj): # need to find proper nearest neighbor
            tstring = "%5s" % (bval+1) 
            newstring = newstring + tstring
#    print alen, mylist
#    print i-j+1,newstring
#    raw_input()
    a2string.append(newstring)
    a2t.append(atype[i])
for i in range (len(pair1)):
    j=pair1[i]
    k=pair2[i]
    tstring = "%5s" % (k+1)
    a2string[j] += tstring
    tstring = "%5s" % (j+1)
    a2string[k] += tstring
    
l2_num=len(a2x)
write_txyz('ctest3.txyz',l2_num,a2t,a2string,a2x,a2y,a2z,ab2t,uc,header)
write_tinker_key('ctest3.key',uc)
print 'Done'