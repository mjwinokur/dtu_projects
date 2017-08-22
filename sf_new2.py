#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  12 15:06:44 2017
A faster way to calculate the structure function by unrolling loops
@author: winokur
"""
import numpy as np
#import math
import time
from mw_library import read_txyz
from mw_library import abc_to_rlp
from mw_library import get_hkl_list_from_cisfile
from mw_library import generate_hkl
from mw_library import sequence_atom_ff
from mw_library import iofiles
from mw_library import write_txyz
from mw_library import write_tinker_key
from rbuildTx import rebuild_Tx
import sys
#
def CH_adjust(aax,aay,aaz,atype,astring,CHdist):
    for i in range(len(atype)):
        if (atype[i] == 'H'):
            ar1 = np.array([aax[i],aay[i],aaz[i]])
            j = int(astring[i])-1 # With hydrogen only one term and the indexing starts at 1 not 0
            ar2 = np.array([aax[j],aay[j],aaz[j]])
            dr = ar1-ar2
            dist = np.sqrt(np.dot(dr,dr))
            ar1 = ar2 +dr*CHdist/dist
            aax[i] = ar1[0]
            aay[i] = ar1[1]
            aaz[i] = ar1[2]
    return aax,aay,aaz
#
#
tfile ='';hkl_output_file='';text=''
tfile,hkl_output_file,text = iofiles(sys.argv[1:])
#print tfile; print hkl_output_file; print text
adjustCH = "F"
keepall = 'false'
tilta = 0.0; tiltb = 0.0; tiltc = 0.0
alpha=0.; beta=0.; gamma = 0.0
efile = 'none'
if (text != ""):
#    print "with special:",text
    alist = text.split()
#    print alist
    j=0
    for i in alist:
#        print j,i
        if (i == 'CHadjust'):
            adjustCH = "T"
            CHdist = float(alist[j+1])
        if (alist[j] == 'tilta'): # An arbitrary space
            tilta = float(alist[j+1])
            print 'tilta: ',tilta
        elif (alist[j] == 'tiltb'): # An arbitrary space
            tiltb = float(alist[j+1])
            print 'tiltb: ',tiltb
        elif (alist[j] == 'tiltc'): # An arbitrary space
            tiltc = float(alist[j+1])
            print 'tiltc: ',tiltc
        elif (alist[j] == 'alpha'): # An arbitrary space
            alpha = float(alist[j+1])
            print 'alpha: ',alpha
        elif (alist[j] == 'beta'): # An arbitrary space
            beta = float(alist[j+1])
            print 'beta: ',beta
        elif (alist[j] == 'gamma'): # An arbitrary space
            gamma = float(alist[j+1])
            print 'gamma: ',gamma
        elif (alist[j] == 'monomer'): # An arbitrary space
            monomer = alist[j+1]
            mono_ct = int(alist[j+2])  # Only use this many monomers in the SF calculation
            print 'monomer: ',monomer  # T2 or T3
        elif (alist[j] == 'export'): # An arbitrary space
            efile = alist[j+1]
            print 'Export: ',efile,'.txyz and ',efile,'.key'
        elif (alist[j] == 'keepall'): # An arbitrary space
            keepall = 'true'
            print 'Keeping all Miller indices: ',
        j += 1

#print 'text:',text
#print 'Argument List:', tfile, hkl_output_file
#raw_input()
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
#Gmax=1.2 # maximum wavevector
dmin=2.*np.pi/Gmax
#tfile="/home/winokur/dtu_projects/NaT2mw.txyz"
#tfile="/home/winokur/Documents/gixrd/NaT2mw.txyz"
hkl_input_file = ''
#tfile="/home/winokur/dtu_projects/ztest.100" # annealing at 300 K
#hkl_output_file = 'ztest.hkl'
if (tfile == '' or hkl_output_file == ''):
    if (tfile == '' and hkl_output_file == ''):
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
        tfile="/home/winokur/dtu_projects/z2test.200" # at 900 K
        hkl_output_file = 'z2test.hkl_200'
        tfile="/home/winokur/dtu_projects/z2test.xyz_3" # at 0 K
        hkl_output_file = 'z2test.hkl_xyz_3'
        tfile="/home/winokur/dtu_projects/z2test.xyz_5" # at 0 K new unit cell
        hkl_output_file = 'z2test.hkl_xyz_5'
    else:
        print('Both the input and the output files are needed!')
        raw_input()
#tfile="/home/winokur/dtu_projects/ztest.xyz_2"
#hkl_output_file = 'ztest2.hkl'
#tfile="/home/winokur/dtu_projects/ztest.xyz"
#hkl_output_file = 'ztest0.hkl'
#tfile="/home/winokur/dtu_projects/wtest.xyz_14"
#hkl_output_file = 'wtest14.hkl'
l2_num, aax,aay,aaz,atype, astring, abt, uc, header = read_txyz(tfile)
if (monomer == 'T2' or monomer == 'T3'):
    if (monomer == 'T2'):
        mnum =96
    elif (monomer == 'T3'):
        mnum =110
    l2_num =  mnum*mono_ct
    print 'Elimating atoms beyond number:',l2_num
    aax = aax[0:l2_num]
    aay = aay[0:l2_num]
    aaz = aaz[0:l2_num]
    del atype[l2_num:]
    del astring[l2_num:]
    del abt[l2_num:]
 
if (l2_num != len(aax)):
    print ' Length mismatch', l2_num, len(aax),' Nominal monomer count',l2_num/mnum
    raw_input()
if (len(uc)== 0):
    raw_input('Stop, this tinker file requires unit cell information')
if (adjustCH == 'T'):
    print " Adjusting CH distance to ",CHdist,"Angstroms"
    aax,aay,aaz = CH_adjust(aax,aay,aaz,atype,astring,CHdist) # 1.10 to 1.00 Angstroms
if (alpha != 0.):
    uc[3]=alpha
if (beta != 0.):
    uc[4]=beta
if (gamma != 0.):
    uc[5]=gamma
if (tilta != 0.0 or tiltb != 0.0 or tiltc!=0):
    print 'Rebuilding'
    if (monomer != 'T2' and monomer != 'T3'):
        print 'Not a valid monomer type for a rebuild:',monomer 
        raw_input()
    aax,aay,aaz = rebuild_Tx(tilta,tiltb,tiltc,uc,aax,aay,aaz,mnum)
if (efile != 'none'):
    write_txyz(efile+'.txyz',l2_num,atype,astring,aax,aay,aaz,abt,uc,header)
    write_tinker_key(efile+'.key',uc)
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
# Inorder to make the python code faster this SF calculation is nearly incomprehensible
#
Rffjk=0.0
for ii in range(iat):  # q = 0 calculation
    [ia0,ia1,ia2,ia3,ia4,ia5,ia6,ia7,ia8]=aff[ii]
    afact[ii] =   ia8 + ia0 +ia2 +ia4 +ia6
for ih in range(l2_num):
    aform[ih] -= 1 # Make the indexing start at zero instead of 1
for ih in range(l2_num):
    ff=afact[aform[ih]]
    Rffjk=Rffjk+ff
norm0 = 100000./(Rffjk*Rffjk) 
Ga = np.zeros(l_num, dtype=float)
Gb = np.zeros(l_num, dtype=float)
Gc = np.zeros(l_num, dtype=float)
for j in range(l_num):
    i=sort_index[m_num-j]  # sort at the beginning
    [Ga[j],Gb[j],Gc[j]]=Ghkl[i]
######################################################
Gaa = np.array(Ga[:],dtype=float)
Gbb = np.array(Gb[:],dtype=float)
Gcc = np.array(Gc[:],dtype=float)
Ax = np.array(aax[:],dtype=float)
Ay = np.array(aay[:],dtype=float)
Az = np.array(aaz[:],dtype=float)
Rphase = np.outer(Gaa,Ax)
Rphase = np.add(Rphase,np.outer(Gbb,Ay))
Rphase = np.add(Rphase,np.outer(Gcc,Az))
Iphase = Rphase
Rphase = np.cos(Rphase)
Iphase = np.sin(Iphase)
# Calculate the atomic form factors
#    G2=np.dot(Ghkl[i],Ghkl[i])
#    sval = -1.*G2*inv4pi2
G2=(np.square(Gaa)+np.square(Gbb)+np.square(Gcc))*(-inv4pi2)
ff = np.empty([iat,l_num])
for ii in range(iat):
    [ia0,ia1,ia2,ia3,ia4,ia5,ia6,ia7,ia8]=aff[ii]
    ff[ii,:] = np.array(ia8 + ia0*np.exp(ia1*G2)+ia2*np.exp(ia3*G2)+ia4*np.exp(ia5*G2)+ia6*np.exp(ia7*G2))
fff = np.empty([l2_num])
# This can probably be automated
if (iat > 12):
    print "iat too big, increase the size of flist in the py file"
    raw_input()
flist = [[],[],[],[],[],[],[],[],[],[],[],[],[]]
j = 0
for i in aform:
    flist[i].append(j)  # Organizes the atom types into an array list
    j += 1
# Above is a little better but there are probably even faster methods
print 'Now for the slow part, have patience'
then = time.time()
#
for ih in range(l_num):
    for j in range(iat):
        temp = ff[j,ih]
        for i in flist[j]:
            fff[i]=temp
# The two line method below is three times slower!
#    for j in range(l2_num):
#        fff[j] = ff[aform[j],ih]
    Rphase[ih,:] = np.multiply(fff,Rphase[ih,:])
    Iphase[ih,:] = np.multiply(fff,Iphase[ih,:])
RSphase = np.sum(Rphase,axis=1)
ISphase = np.sum(Iphase,axis=1)
#    Iq_temp = (Rffjk*Rffjk+Iffjk*Iffjk)*norm0
RSphase = np.square(RSphase)
ISphase = np.square(ISphase)
RSphase = np.add(RSphase,ISphase)*norm0
n_keep = 0; n_throw = 0; Iq = []; Iq_idx = []
for j in range(l_num):
    i = sort_index[m_num-j]
    Iq_temp = RSphase[j]
    if (Iq_temp < 0.001 and keepall == 'false'): 
        n_throw += 1  # Don't keep essentially zero intensities
    else:
        Iq.append(Iq_temp)
        Iq_idx.append(i)
        n_keep += 1
#    if (500*int(j/500) == j):
#        print '\n'
now = time.time()
print "Slow part elapsed time ",now - then,' seconds'
print ' SF done with ',j,'/',l_num,' of which ',n_keep,' reflections were kept'
#
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
# These rings are a 2D powder average
#        renorm2 = np.sqrt(Gpara*Gpara+Gperp)
        renorm2 = Gpara
        if (Gpara < 0.05): #  Along the perp
#            renorm2=np.sqrt(Gperp)
            renorm2=1.
        mynewstring="%6s%7s%7s%13.6f%14.4f%5.0f" %  (hv[i], kv[i], lv[i], d[i], Iq[j]*renorm/renorm2, mult[i] ) # Ignore (000) reflection
    else:
        mynewstring="%6s%7s%7s%13.6f%14.4f%5.0f" %  (hv[i], kv[i], lv[i], d[i], Iq[j]*renorm, mult[i]) 
    f.write(mynewstring+'\n')

f.close()  # close write file
print 'sf_new2 fini'