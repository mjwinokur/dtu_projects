#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 22:09:16 2017

@author: winokur
"""
# This file takes a list of tinker txyz files and averages the structure factors
import numpy as np
#import math
#import time
from mw_library import read_txyz_info_uc
from mw_library import iofiles
from mw_library import read_hklI_full_file
import sys
import os
txyz_name_base = '/home/winokur/dtu_projects/S1A_101_15n/NaT2101a'
i_start=81
i_stop=160
home='/home/winokur/dtu_projects/'
text = 'CHadjust 1.0 tiltb 0. tilta 0 tiltc 0 monomer T2 72 keepall silent'  
text2 = '"CHadjust 1.0 tiltb 0. tilta 0 tiltc 0 monomer T2 72 keepall silent"'  
txyz_name_base,hkl_output_file,text = iofiles(sys.argv[1:])
tinkfile = 'temp.txyz'
hklfile = 'temp.hkl'
alpha=0.;beta=0.;gamma=0.
#Example python /home/winokur/dtu_projects/smush_hkl.py -i /home/winokur/dtu_projects/S1A_101_15n/NaT2101a -o temp_ave.hkl -t "start 81 stop 82 monomer T2 count 72"
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
        elif (alist[j] == 'stop'): 
            i_stop = int(alist[j+1])
            print 'stop: ',i_stop
        elif (alist[j] == 'monomer'): 
            monomer = alist[j+1]
            print 'monomer: ',monomer
        elif (alist[j] == 'count'): 
            mono_ct = int(alist[j+1])
            print 'count: ',mono_ct
        elif (alist[j] == 'alpha'): 
            alpha = float(alist[j+1])
            print 'alpha: ', alpha
        elif (alist[j] == 'beta'): 
            beta = float(alist[j+1])
            print 'beta: ', beta
        elif (alist[j] == 'gamma'): 
            gamma = float(alist[j+1])
            print 'gamma: ', gamma
        j += 1
#
norm=1./float(i_stop-i_start+1)
str_num=str(i_start)
str_len=len(str_num)
if (str_len == 1):
    tfile=txyz_name_base+'.00'+str_num
elif (str_len == 2):
    tfile=txyz_name_base+'.0'+str_num
elif (str_len == 3):
    tfile=txyz_name_base+'.'+str_num
s1_cmd= 'cp '+tfile+' temp.txyz'
print s1_cmd
os.system(s1_cmd)
uc = read_txyz_info_uc(tinkfile)
if (len(uc)== 0):
    raw_input('Stop, this tinker file requires unit cell information')
if (alpha != 0.):
    uc[3]=alpha
if (beta != 0.):
    uc[4]=beta
if (gamma != 0.):
    uc[5]=gamma
# -i NaT2101a_81_160.txyz -o NaT2101a_81_160.hkl -t 
s2_cmd= 'python '+home+'sf_new2.py -i '+tinkfile+' -o '+hklfile+' -t '+text2
print s2_cmd
os.system(s2_cmd)
# This version does NOT multiply intensity*multiplicity!!!!
hkl,I0,dspace,mult = read_hklI_full_file(hklfile)
I0=np.array(I0)
l_num = len(hkl)
#
for i in range(i_start+1,i_stop+1):
    str_num=str(i)
    str_len=len(str_num)
    if (str_len == 1):
        tfile=txyz_name_base+'.00'+str_num
    elif (str_len == 2):
        tfile=txyz_name_base+'.0'+str_num
    elif (str_len == 3):
        tfile=txyz_name_base+'.'+str_num
    s1_cmd= 'cp '+tfile+' temp.txyz'
    print '1 Executing :',s1_cmd
    os.system(s1_cmd)
    print '2 Executing :',s2_cmd
    os.system(s2_cmd)
#    raw_input('Paused')
# -i NaT2101a_81_160.txyz -o NaT2101a_81_160.hkl -t 
# Now open the hkl file
# This version does NOT multiply intensity*multiplicity!!!!
    hkl,I,dspace,mult = read_hklI_full_file(hklfile)
    l_num_new = len(hkl)
    if (l_num_new != l_num):
        print(' The number of Miller indicies has changed, this should not occur')
        raw_input()
    else:
        l_num = l_num_new
    I=np.array(I)
    I0 += I
Iq = I0*norm
# Now to save the final result
f = open(hkl_output_file,'w')
f.write('     h      k      l  d-spacing       F^2     multiplicity'+'\n')
#for j in range(1,l_num): # exclude 0 0 0 reflection
for i in range(len(Iq)): # exclude 0 0 0 reflection
#    print 'hkl',i,hkl[i],dspace[i],Iq[i]
    [h,k,l] = hkl[i]
    mynewstring="%6s%7s%7s%13.6f%14.4f%5.0f" %  (int(h),int(k),int(l), dspace[i], Iq[i], mult[i]) 
    f.write(mynewstring+'\n')
f.close()  # close write file
print 'smush_hkl fini'