#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 11:12:26 2017

@author: winokur
"""
# smush_txyz.py is designed to open a series of MD simulations and average the atom positon
import numpy as np
#import math
#import time
from mw_library import read_txyz
#from mw_library import iofiles
from mw_library import write_txyz
from mw_library import write_tinker_key
from mw_library import iofiles
import sys
#
def find_center(aax,aay,aaz,monomer):
    datamean=[]
    if (monomer == 'T2'):
        xt = np.append(np.asarray(aax[0:24],dtype=float),np.asarray(aax[48:72],dtype=float))
        yt = np.append(np.asarray(aay[0:24],dtype=float),np.asarray(aay[48:72],dtype=float))
        zt = np.append(np.asarray(aaz[0:24],dtype=float),np.asarray(aaz[48:72],dtype=float))
        data = np.concatenate((xt[:,np.newaxis],yt[:,np.newaxis],zt[:,np.newaxis]),axis=1)
# Calculate the mean of the points, i.e. the 'center' of the cloud
        datamean.append(data.mean(axis=0))
# Do an SVD on the mean-centered data.
#
        xt = np.append(np.asarray(aax[24:48],dtype=float),np.asarray(aax[72:96],dtype=float))
        yt = np.append(np.asarray(aay[24:48],dtype=float),np.asarray(aay[72:96],dtype=float))
        zt = np.append(np.asarray(aaz[24:48],dtype=float),np.asarray(aaz[72:96],dtype=float))
        data = np.concatenate((xt[:,np.newaxis],yt[:,np.newaxis],zt[:,np.newaxis]),axis=1)
# Calculate the mean of the points, i.e. the 'center' of the cloud
        datamean.append(data.mean(axis=0))
    elif (monomer == 'T3'):
        xt = np.asarray(aax[0:55],dtype=float)
        yt = np.asarray(aay[0:55],dtype=float)
        zt = np.asarray(aaz[0:55],dtype=float)
        data = np.concatenate((xt[:,np.newaxis],yt[:,np.newaxis],zt[:,np.newaxis]),axis=1)
#
        xt = np.asarray(aax[55:110],dtype=float)
        yt = np.asarray(aay[55:110],dtype=float)
        zt = np.asarray(aaz[55:110],dtype=float)
        data = np.concatenate((xt[:,np.newaxis],yt[:,np.newaxis],zt[:,np.newaxis]),axis=1)
# Calculate the mean of the points, i.e. the 'center' of the cloud
        datamean.append(data.mean(axis=0))
    return datamean
#
def rewrap(monomer):
    if (monomer == 'T2'):
        aone=np.ones(24)
        azero=np.zeros(24)
        a0 = aone
        a0 = np.append(a0,azero)
        a0 = np.append(a0,aone)
        a0 = np.append(a0,azero)
        a1 = azero
        a1 = np.append(a1,aone)
        a1 = np.append(a1,azero)
        a1 = np.append(a1,aone)
    elif (monomer == 'T3'):
        aone=np.ones(55)
        azero=np.zeros(55)
        a0 = aone
        a0 = np.append(a0,azero)
        a1 = azero
        a1 = np.append(a1,aone)
    return a0,a1
#
# monomer T2 72
mono_ct = 72
monomer = 'T2'
txyz_name_base = '/home/winokur/dtu_projects/S1A_101_15n/NaT2101a'
i_start=81
i_stop=160
txyz_name_base,hkl_output_file,text = iofiles(sys.argv[1:])
#Example python /home/winokur/dtu_projects/smush_txyz -i /home/winokur/dtu_projects/S1A_101_15n/NaT2101a -t "start 81 stop 160 monomer T2 count 72"
if (text != ""):
#    print "with special:",text
    alist = text.split()
#    print alist
    j=0
    for i in alist:
#        print j,i
        if (alist[j] == 'start'):
            i_start = int(alist[j+1])
            print 'start: ',i_start
        elif (alist[j] == 'stop'): # An arbitrary space
            i_stop = int(alist[j+1])
            print 'stop: ',i_stop
        elif (alist[j] == 'monomer'): # An arbitrary space
            monomer = alist[j+1]
            print 'monomer: ',monomer
        elif (alist[j] == 'count'): # An arbitrary space
            mono_ct = int(alist[j+1])
            print 'count: ', mono_ct
        j += 1
ict = 1
if (monomer == 'T2' or monomer == 'T3'):
    if (monomer == 'T2'):
        mnum =96
    elif (monomer == 'T3'):
        mnum =110
str_num=str(i_start)
str_len=len(str_num)
if (str_len == 1):
    tfile=txyz_name_base+'.00'+str(i_start)
elif (str_len == 2):
    tfile=txyz_name_base+'.0'+str(i_start)
elif (str_len == 3):
    tfile=txyz_name_base+'.'+str(i_start)
l2_num, aax0,aay0,aaz0,atype, astring, abt, uc, header = read_txyz(tfile)
l2_num =  mnum*mono_ct
print 'Elimating atoms beyond number:',l2_num
taax = np.array(aax0[0:l2_num])
taay = np.array(aay0[0:l2_num])
taaz = np.array(aaz0[0:l2_num])
del atype[l2_num:]
del astring[l2_num:]
del abt[l2_num:]
datamean0=[]
for i in range(mono_ct):
    mi=mnum*i
    mf=mi+mnum
    datamean0.append(find_center(taax[mi:mf],taay[mi:mf],taaz[mi:mf],monomer))
#
a0,a1 = rewrap(monomer)
#
for i in range(i_start+1,i_stop+1):
    str_num=str(i)
    str_len=len(str_num)
    if (str_len == 1):
        tfile=txyz_name_base+'.00'+str(i)
    elif (str_len == 2):
        tfile=txyz_name_base+'.0'+str(i)
    elif (str_len == 3):
        tfile=txyz_name_base+'.'+str(i)
    l2_num2, aaxx,aayy,aazz,atype, astring, abt, uc, header = read_txyz(tfile)
    aax = np.array(aaxx[0:l2_num])
    aay = np.array(aayy[0:l2_num])
    aaz = np.array(aazz[0:l2_num])
# Positions can be off by one unit cell
    ucx=uc[0]*0.5
    ucy=uc[1]*0.5
    ucz=uc[2]*0.5
    datamean=[]
    for j in range(mono_ct):
        mi=mnum*j
        mf=mi+mnum
        datamean.append(find_center(aax[mi:mf],aay[mi:mf],aaz[mi:mf],monomer))
    datadiff = np.subtract(datamean0,datamean)
#    datadiff = np.subtract(datamean,datamean0)
# Check and adjust if a rewrap is needed
    for j in range(mono_ct):
        mi=mnum*j
        mf=mi+mnum
#
        x0 = np.abs(datadiff[j,0,0])
        if (x0 > ucx):
            amult = np.float(np.sign(datadiff[j,0,0]))
            temp = amult*uc[0]*a0
            aax[mi:mf] = np.add(aax[mi:mf],temp)
#            print 'Rewrap x0 for i,j:',i,j
#            
        y0 = np.abs(datadiff[j,0,1])
        if (y0 > ucy):
            amult = np.float(np.sign(datadiff[j,0,1]))
            temp = amult*uc[1]*a0
            aay[mi:mf] = np.add(aay[mi:mf],temp)
#            print 'Rewrap y0 for i,j:',i,j
#
        z0 = np.abs(datadiff[j,0,2])
        if (z0 > ucz):
            amult = np.float(np.sign(datadiff[j,0,2]))
            temp = amult*uc[2]*a0
            aaz[mi:mf] = np.add(aaz[mi:mf],temp)
#            print 'Rewrap z0 for i,j:',i,j
#
        x1 = np.abs(datadiff[j,1,0])
        if (x1 > ucx):
            amult = np.float(np.sign(datadiff[j,1,0]))
            temp = amult*uc[0]*a1
            aax[mi:mf] = np.add(aax[mi:mf],temp)
#            print 'Rewrap x1 for i,j:',i,j
#
        y1 = np.abs(datadiff[j,1,1])
        if (y1 > ucy):
            amult = np.float(np.sign(datadiff[j,1,1]))
            temp = amult*uc[1]*a1
            aay[mi:mf] = np.add(aay[mi:mf],temp)
#            print 'Rewrap y1 for i,j:',i,j
#
        z1 = np.abs(datadiff[j,1,2])
        if (z1 > ucz):
            amult = np.float(np.sign(datadiff[j,1,2]))
            temp = amult*uc[2]*a1
#            print 'Base', j,mi,mf,uc[2]
#            print aaz0[mi:mf]
#            print 'Before'
#            print aaz[mi:mf]
            aaz[mi:mf] = np.add(aaz[mi:mf],temp)
#            print 'After'
#            print aaz[mi:mf]
#            print 'Rewrap z1 for i,j:',i,j
#            raw_input()
#    datamean=[]
#    for j in range(mono_ct):
#        mi=mnum*j
#        mf=mi+mnum
#        datamean.append(find_center(aax[mi:mf],aay[mi:mf],aaz[mi:mf],monomer))
#    datadiff = np.subtract(datamean,datamean0)
#    print i
#    raw_input()
    taax += aax
    taay += aay
    taaz += aaz
    ict += 1
#    raw_input()
a =1./float(ict)
aax = taax * a
aay = taay * a
aaz = taaz * a
txyz_out=txyz_name_base+'_'+str(i_start)+'_'+str(i_stop)
write_txyz(txyz_out+'.txyz',l2_num,atype,astring,aax,aay,aaz,abt,uc,header)
write_tinker_key(txyz_out+'.key',uc)
print 'Fini'
# end