#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 09:59:03 2017

@author: winokur
"""
# This function calculates the average C-C and C-S pair distances
from __future__ import print_function
import numpy as np
#import math
#import time
from mw_library import read_txyz
#from mw_library import abc_to_rlp
#from mw_library import get_hkl_list_from_cisfile
#from mw_library import generate_hkl
#from mw_library import sequence_atom_ff
from mw_library import iofiles
#from mw_library import write_txyz
#from mw_library import write_tinker_key
#from rbuildTx import rebuild_Tx
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
mpl.rcParams['legend.fontsize'] = 10
#import matplotlib.gridspec as gridspec
##from matplotlib import rc
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
def pairs_from_tinker(atype,astring,aax,aay,aaz):
    pairs = [['C','C'],['C','S']]
    nbr_idx = []
    for i in range(len(aax)):
        line = astring[i]
        newline = ' '.join(line.split())
        mylist = newline.split(" ")
#        print(mylist[0])
        idx = int(mylist[0])-1
        atm1=atype[idx]
        atm0=atype[i]
        for j in pairs:
#            print(j)
            if (atm0 == j[0] and atm1 == j[1]):
                if (idx < i):
                    nbr_idx.append([idx,i])
                else:
                    nbr_idx.append([i,idx])
            elif (atm0 == j[1] and atm1 == j[0]):
                if (idx < i):
                    nbr_idx.append([idx,i])
                else:
                    nbr_idx.append([i,idx])
#        print(atm0,atm1)
#        print(idx,i)
    nbr_idx.sort()  
    lst =[]
    tmp = [0,0]
# This oversamples and so double counting must be eliminated
    for i in nbr_idx:
        if (i != tmp):
            lst.append(i)
            tmp = i
#    print(lst)
    nbr_idx = lst
    xyz = np.asanyarray(3)
    rad = []
#    print(nbr_idx)
# Now to generate pair distances
    for i in nbr_idx:
        [ia,ib] = i
        xyz = [aax[ia]-aax[ib],aay[ia]-aay[ib],aaz[ia]-aaz[ib]]
        rad.append(np.sqrt(np.dot(xyz,xyz)))
    return nbr_idx,rad
#
def pairs_from_radius(atype,aax,aay,aaz,radius):
    atom_types = ['C','S']
    mseq = []
    for i in range(len(aax)):
        for j in atom_types:
            if (atype[i]==j):
                mseq.append(i)
    nbr_idx = []
    xyz = np.asanyarray(3)
    ict = len(mseq)
    for i in range(ict):
        for j in range(i+1,ict):
            ia = mseq[i]
            ib = mseq[j]
            xyz = [aax[ia]-aax[ib],aay[ia]-aay[ib],aaz[ia]-aaz[ib]]
            rad0 = np.sqrt(np.dot(xyz,xyz))
            if (rad0 < radius):
                nbr_idx.append([ia,ib,rad0])
    nbr_idx.sort()  
    lst =[]
    rad = []
    tmp = [0,0,0]
# This oversamples and so double counting must be eliminated
    for i in nbr_idx:
        if (i != tmp):
            lst.append([i[0],i[1]])
            rad.append(i[2])
            tmp = i
#    print(lst)
    nbr_idx = lst
#    print(nbr_idx)
# Now to generate pair distances
    return nbr_idx,rad

input_file ='NaT2mw.txyz'
input_file ='NaT3_single_xtal_all_flip_gap/NaT3_833_allflip_gap_css2.040'
output_file='';text=''
monomer='T2'
monomer='T3'
# Example use 
# python pairdist.py -i NaT3_single_xtal_all_flip_gap/NaT3_833_allflip_gap_css2.040 -t 'mono T3'
ifile,ofile,text = iofiles(sys.argv[1:])
if (ifile != ''):
    input_file =ifile
#print(tfile; print(hkl_output_file; print(text
if (text != ""):
#    print("with special:",text
    alist = text.split()
#    print(alist
    j=0
    for i in alist:
#        print(j,i
        if (alist[j] == 'tilta'): # An arbitrary space
            tilta = float(alist[j+1])
            print('tilta: ',tilta)
        elif (alist[j] == 'mono'): # An arbitrary space
            monomer = alist[j+1]
            print('Monomer type: ',monomer)
        j += 1
#print('text:',text
#print('Argument List:', tfile, hkl_output_file
#raw_input()
l2_num, aax,aay,aaz,atype, astring, abt, uc, header = read_txyz(input_file)
if (monomer == 'T2' or monomer == 'T3'):
    if (monomer == 'T2'):
        mnum =96
    elif (monomer == 'T3'):
        mnum =110
#
#nbr_idx,rad = pairs_from_tinker(atype,astring,aax,aay,aaz) # but this algorithm doesn't get all the pairs
nbr_idx,rad = pairs_from_radius(atype,aax[0:mnum],aay[0:mnum],aaz[0:mnum],2.0) # but this algorithm doesn't get all the pairs

#i1 =[]; i2 = []; npairs=len(rad)
npairs=len(rad)
i1 =np.zeros(npairs,dtype=int); i2 = np.zeros(npairs,dtype=int) 
radii = np.zeros(npairs)
radsqr = np.zeros(npairs)
xyz = np.zeros(3)
t_mean = [];x_m=np.zeros(npairs);y_m=np.zeros(npairs);z_m=np.zeros(npairs);s_t=[]
for j in range(npairs):
    [ia,ib] = nbr_idx[j]
    i1[j]=ia;i2[j]=ib
#    i1.append(ia)
#    i2.append(ib)
    radii[j]=rad[j]
    radsqr[j]=rad[j]*rad[j]
    x_m[j]=(0.5*(aax[ia]+aax[ib]))
    y_m[j]=(0.5*(aay[ia]+aay[ib]))
    z_m[j]=(0.5*(aaz[ia]+aaz[ib]))
#        
n_mono = int(l2_num/mnum)
print('Number of cells is:',n_mono)
#n_mono = 30
for j in range(1,n_mono):
    ic = mnum*j
    i11 = i1+ic
    i22 = i2+ic
    for k in range(npairs):
        xyz = [aax[i11[k]]-aax[i22[k]],aay[i11[k]]-aay[i22[k]],aaz[i11[k]]-aaz[i22[k]]]
        tmp = np.dot(xyz,xyz)
#        print(j,k,i1[k],i11[k],i2[k],i22[k],rad[k],tmp)
#        raw_input()
        radsqr[k]+= tmp
        radii[k] += np.sqrt(tmp)
if (n_mono > 1):
    amult=1./float(n_mono)
    for j in range(npairs):
        tmp = np.asarray([radii[j]*amult])
        tmp2 = np.asarray([radsqr[j]*amult])
        std = np.sqrt(tmp2-tmp*tmp)
#        print(j,'std',tmp,tmp2,std)
#        raw_input()
        s_text = "%5.3f" % (tmp)
        s_text2 = "%4.3f" % (std)
        s_t.append(s_text+'('+s_text2+')')
else:    
    for j in range(npairs):
        s_text = "%5.3f" % (radii[j])
        s_t.append(s_text)
#    print(j,x_m[j],y_m[j],z_m[j],' ',s_t[j])
#
# Now to calculate avg CS and CC
radCC=0.
iCC=0.
radCS=0.
iCS=0.
j=0
for i in nbr_idx:
    if (atype[i[0]] == 'C' and atype[i[1]] == 'C'):
        radCC += radii[j]
        iCC += 1.
    elif ((atype[i[0]] == 'S' and atype[i[1]] == 'C') or (atype[i[0]] == 'C' and atype[i[1]] == 'S')):
        radCS += radii[j]
        iCS += 1.
    j += 1
print('Avg. CC bond:',amult*radCC/iCC,'  Avg. CS bond: ',amult*radCS/iCS)
#
fig = plt.figure()
ax = fig.gca(projection='3d')
#
datasets = [{"x":[aax[i1[i]],aax[i2[i]]], "y":[aay[i1[i]],aay[i2[i]]], "z":[aaz[i1[i]],aaz[i2[i]]], "colour": "black"} for i in range(len(nbr_idx))]

#Axes3D.text(x, y, z, s, zdir=None, **kwargs)
for dataset in datasets:
#    print('here',dataset["x"], dataset["y"], dataset["z"])
    ax.plot(dataset["x"], dataset["y"], dataset["z"], color=dataset["colour"])
for i in range(npairs):   
    ax.text(x_m[i],y_m[i],z_m[i],s_t[i])
ax.legend()
plt.show()
#raw_input()
print('Fini')        