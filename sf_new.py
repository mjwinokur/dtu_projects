#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 15:06:44 2017

@author: winokur
"""
import numpy as np
import math
#import time
from mw_library import read_txyz
from mw_library import abc_to_rlp
from mw_library import get_hkl_list_from_cisfile
from mw_library import generate_hkl
from mw_library import sequence_atom_ff
import sys
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)
raw_input()
##################################################################################
#
# A basic structure function python module
#
###################################################################################
#
# Inputs  
# Gmax maximum wave vector
# Tinker xyz file with positions and unit
# List of hkl indices
#
# Outputs
# hkl file
#
Gmax=3.2 # maximum wavevector
dmin=2.*np.pi/Gmax
#tfile="/home/winokur/dtu_projects/NaT2mw.txyz"
#tfile="/home/winokur/Documents/gixrd/NaT2mw.txyz"
hkl_input_file = ''
#tfile="/home/winokur/dtu_projects/ztest.100" # annealing at 300 K
#hkl_output_file = 'ztest.hkl'
tfile="/home/winokur/dtu_projects/ztest.300" # at 400 K
hkl_output_file = 'ztest.hkl_300'
tfile="/home/winokur/dtu_projects/ztest.500" # at 500 K
hkl_output_file = 'ztest.hkl_500'
tfile="/home/winokur/dtu_projects/ztest.700" # at 550 K
hkl_output_file = 'ztest.hkl_700'
tfile="/home/winokur/dtu_projects/ztest.900" # at 600 K
hkl_output_file = 'ztest.hkl_900'
tfile="/home/winokur/dtu_projects/ztest.1100" # at 700 K
hkl_output_file = 'ztest.hkl_1100'
tfile="/home/winokur/dtu_projects/ztest.1300" # at 800 K
hkl_output_file = 'ztest.hkl_1300'
tfile="/home/winokur/dtu_projects/ztest.1500" # at 900 K
hkl_output_file = 'ztest.hkl_1500'
tfile="/home/winokur/dtu_projects/z2test.030" # at 900 K
hkl_output_file = 'z2test.hkl_30'
#tfile="/home/winokur/dtu_projects/ztest.xyz_2"
#hkl_output_file = 'ztest2.hkl'
#tfile="/home/winokur/dtu_projects/ztest.xyz"
#hkl_output_file = 'ztest0.hkl'
#tfile="/home/winokur/dtu_projects/wtest.xyz_14"
#hkl_output_file = 'wtest14.hkl'
l2_num, aax,aay,aaz,atype, astring, abt, uc, header = read_txyz(tfile)
if (len(uc)== 0):
    raw_input('Stop, this tinker file requires unit cell information')
# From a,b,c, alpha, beta,gamma get the x,y,z components
astar,bstar,cstar = abc_to_rlp(uc)
if (hkl_input_file != ''):
    hv,kv,lv,d,Ghkl,mult,sort_index = get_hkl_list_from_cisfile(hkl_input_file,astar,bstar,cstar,Gmax)
else:
    hv,kv,lv,d,Ghkl,mult,sort_index = generate_hkl(Gmax,astar,bstar,cstar)
l_num = len(hv)
astarhat=astar/np.sqrt(np.dot(astar,astar))
#
# Find out which atoms are being used
#
aform,aff = sequence_atom_ff(atype)
iat = len(aff)
afact = np.zeros(iat) # Define an empty array
#
inv4pi2=1./(16.*np.pi*np.pi)
m_num = l_num - 1
#
# Normalize the I(000) reflection to 100000.
#
Rffjk=0.0
for ii in range(iat):
    [ia0,ia1,ia2,ia3,ia4,ia5,ia6,ia7,ia8]=aff[ii]
    afact[ii] =   ia8 + ia0 +ia2 +ia4 +ia6
for ih in range(l2_num):
#    Rj = np.array([aax[ih],aay[ih],aaz[ih]]) # doesn't save time
    ff=afact[aform[ih]-1]  # Indexing starts at one
    Rffjk=Rffjk+ff
norm0 = 100000./(Rffjk*Rffjk) 
#
Ga = np.zeros(l_num, dtype=float)
Gb = np.zeros(l_num, dtype=float)
Gc = np.zeros(l_num, dtype=float)
for j in range(l_num):
    i=sort_index[m_num-j]
    [Ga[j],Gb[j],Gc[j]]=Ghkl[i]
#
Rj = []
for ih in range(l2_num):
    Rj.append([aax[ih],aay[ih],aaz[ih]])
    aform[ih] -= 1 # Make they indexing start at zero
######################################################
Iq = []
Iq_idx = []
n_keep = 0
n_throw = 0
for j in range(l_num):
    i=sort_index[m_num-j]
    Rffjk=0.0
    Iffjk=0.0
    Gaa =Ga[j];Gbb =Gb[j]; Gcc =Gc[j]
    G2=Gaa*Gaa+Gbb*Gbb+Gcc*Gcc
#    G2=np.dot(Ghkl[i],Ghkl[i])
    sval = -1.*G2*inv4pi2
# Calculate the atomic form factors
    for ii in range(iat):
        [ia0,ia1,ia2,ia3,ia4,ia5,ia6,ia7,ia8]=aff[ii]
        afact[ii] =   ia8 + ia0*math.exp(ia1*sval)+ia2*math.exp(ia3*sval)
        afact[ii]=afact[ii]+ia4*math.exp(ia5*sval)+ia6*math.exp(ia7*sval)
# Calculate the structure function
    for ih in range(l2_num):
#        phase=np.dot(Ghkl[i],Rjj)
#        phase=np.dot(Ghkl[i],[aax[ih],aay[ih],aaz[ih]])
        phase=Gaa*aax[ih]+Gbb*aay[ih]+Gcc*aaz[ih]
        ff=afact[aform[ih]]  # Indexing starts at zero now
        Rffjk += ff * math.cos(phase)  # 30% execution for cos
        Iffjk -= ff * math.sin(phase)  # 30% execution for sin
    Iq_temp = (Rffjk*Rffjk+Iffjk*Iffjk)*norm0
    if (Iq_temp < 0.001 ): 
        n_throw += 1  # Don't keep essentially zero intensities
    else:
        Iq.append(Iq_temp)
        Iq_idx.append(i)
        n_keep += 1
    if (500*int(j/500) == j):
        print '\n'
        print ' SF at ',j,'/',l_num,' of which ',n_keep,' reflections were kept'
#        print(time.time(), time.clock())
# Between numpy and math plus some changes in the array references gave a factor of two speed up
#
#for j in range(l_num):
#    i=sort_index[m_num-j]
#    mynewstring="%5s%5s%5s%13.6f%14.4f%5.0f" %  (hv[i], kv[i], lv[i], d[i], Iq[j],mult[i])
#    print mynewstring
#    print j,i,hv[i],kv[i],lv[i],np.sqrt(np.dot(Ghkl[i],Ghkl[i])),d[i],Iq[j]
#    raw_input()
powder2D = 'true'
renorm = 1.0
f = open(hkl_output_file,'w')
f.write('     h      k      l  d-spacing       F^2     multiplicity'+'\n')
#for j in range(1,l_num): # exclude 0 0 0 reflection
for j in range(1,len(Iq)): # exclude 0 0 0 reflection
    i=Iq_idx[j]    
    if (powder2D == 'true'):
        Gperp = abs(np.dot(Ghkl[i],astarhat))
        Gmag2 = np.dot(Ghkl[i],Ghkl[i])
#        print j,Gmag2,Gperp*Gperp
#        raw_input()
        Gpara= np.sqrt(abs(Gmag2-Gperp*Gperp))        
        renorm2 = np.sqrt(Gpara*Gpara+Gperp)
        if (Gpara < 0.05): #  Along the perp
            renorm2=np.sqrt(Gperp)
        mynewstring="%6s%7s%7s%13.6f%14.4f%5.0f" %  (hv[i], kv[i], lv[i], d[i], Iq[j]*renorm/renorm2, mult[i] ) # Ignore (000) reflection
    else:
        mynewstring="%6s%7s%7s%13.6f%14.4f%5.0f" %  (hv[i], kv[i], lv[i], d[i], Iq[j]*renorm, mult[i]) 
    f.write(mynewstring+'\n')

f.close()  # close write file
print 'fini'