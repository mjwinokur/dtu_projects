#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 1 11:27:47 2017

@author: winokur
"""
#import os
#import openbabel, re  #used in obtest
import numpy as np
#import math
#from gsd.hoomd import Snapshot,open,HOOMDTrajectory
# Without these two cryptic lines below the "import gsd" breaks
#import hoomd
#from hoomd import md
from gsd.hoomd import Snapshot
z = Snapshot()
import gsd.pygsd  # line below only gets hoomd and not pygsd 
import gsd
from mw_library import read_txyz
from mw_library import read_txyz_strip_uc 
from mw_library import tilt_calc 
#from mw_library import write_txyz
from mw_library import write_txyz_2
from mw_library import io_hoomd_params
from mw_library import olig_recenter
from imp_calc_ang import improper_calc_ang
from imp_calc_ang import angle_between
from obtest import btd_info
from forcefield_decode_mm import mm3_decode 
# Decode a tinker file
#mfile = 'NaT2mw.txyz'
#mfile = 'ytest.txyz'
mfile = 'ztest.xyz_2'
temp_file ='temp.txyz'
uc = read_txyz_strip_uc(mfile,temp_file)
if (len(uc)== 0):
    raw_input('Stop, this tinker file requires unit cell information')
l2_num, aax,aay,aaz,atype, astring, abt, uc2, header = read_txyz(temp_file)
if (len(uc2) != 0):
    raw_input('Stop and remome txyz file unit cell line, openbabel interpreter breaks otherwise')
#
#write_txyz('wtest.xyz',l2_num,atype,astring,aax,aay,aaz,abt,uc,header)
#uc = [20.502, 5.989, 8.07, 90., 96.8549, 90.]
# l2_num number of atoms in positions aax,aay,aaz
# atype is the letter element
# astring is the residual bonding information if needed
# abt is currently the MM3 number that identify the bonding configuration
#
# Need unit cell infomation and pair bonding information obabel   
# And force field information
#
fffile ="/home/winokur/MolecularTools/tinker/params/mm3.prm"  # MM force constants
# Now to find what kind of atoms we have
ptypes = [] # link list to coordinate atom/bond type
aatypes = []
for i in abt:
    if i not in ptypes:
        ptypes.append(i)
for i in atype:
    if i not in aatypes:
        aatypes.append(i)
tid = [(ptypes.index(i)) for i in abt] # e.g. encode from 42(S) 2(C) 5(H) to 0 1 2
pair1,pair2,bangle1,bangle2,bangle3,tangle1,tangle2,tangle3,tangle4,mol = btd_info("txyz",temp_file)
pair1 = [(pair1[i]-1) for i in range(len(pair1))]
pair2 = [(pair2[i]-1) for i in range(len(pair2))]
bangle1 = [(bangle1[i]-1) for i in range(len(bangle1))] # indexing needs to start at zero
bangle2 = [(bangle2[i]-1) for i in range(len(bangle2))]
bangle3 = [(bangle3[i]-1) for i in range(len(bangle3))]
tangle1 = [(tangle1[i]-1) for i in range(len(tangle1))]
tangle2 = [(tangle2[i]-1) for i in range(len(tangle2))]
tangle3 = [(tangle3[i]-1) for i in range(len(tangle3))]
tangle4 = [(tangle4[i]-1) for i in range(len(tangle4))]
#print len(pair1),len(bangle1),len(tangle1)
# So we can have S-C or C-S but we use the same force field parameters
bdtypes = []
bdtypeid = []
bdseq = []
angseq = []
torseq = []
for i in range(len(pair1)):
    temp =  atype[pair1[i]]+str(abt[pair1[i]])+atype[pair2[i]]+str(abt[pair2[i]])  # shift to get proper reference
    bdseq.append(temp)
    if temp not in bdtypes:
        bdtypes.append(temp) # all of these are preliminary because the parameterization may have multiple terms
    bdtypeid.append(bdtypes.index(temp))  # all of these are preliminary because the parameterization may have multiple terms  
angtypes = []
angtypeid = []
for i in range(len(bangle1)):
    temp =  atype[bangle1[i]]+str(abt[bangle1[i]])+atype[bangle2[i]]+str(abt[bangle2[i]])+atype[bangle3[i]]+str(abt[bangle3[i]])  # shift to get proper reference
    angseq.append(temp)
    if temp not in angtypes:
        angtypes.append(temp) # all of these are preliminary because the parameterization may have multiple terms
    angtypeid.append(angtypes.index(temp))    # all of these are preliminary because the parameterization may have multiple terms
tortypes = []
tortypeid = []
for i in range(len(tangle1)):
    temp =  atype[tangle1[i]]+str(abt[tangle1[i]])+atype[tangle2[i]]+str(abt[tangle2[i]])+atype[tangle3[i]]+str(abt[tangle3[i]])+atype[tangle4[i]]+str(abt[tangle4[i]])  # shift to get proper reference
    torseq.append(temp)
    if temp not in tortypes:
        tortypes.append(temp)  # all of these are preliminary because the parameterization may have multiple terms
    tortypeid.append(tortypes.index(temp)) #   all of these are preliminary because the parameterization may have multiple terms 
# Now extract the MM3 parameters with the potential bonding of the specific molecular structure here
amass,bdtypes_tot,bond1,bond2,angtypes_tot,angk,angt0,tortypes_tot,tor1,tor2,tor3,vdw1,vdw2,vdwpr1,vdwpr2,dipole_par,dipole_types = mm3_decode(fffile,ptypes,aatypes,bdtypes,angtypes,tortypes)
tor4 = []
for i in range(len(tor1)):
    tor4.append(0.0)
angt0 = map(np.deg2rad, angt0) # convert degrees to radians
#
# Now to deduce molecules
# This assumes that openbabel produces a sequenced numbering of the pairs
# So far, so good
mol_ct=1
mol_idx = []
for i in range(l2_num):
    mol_idx.append(-2)
mol_idx[pair1[0]]=mol_ct
mol_idx[pair2[0]]=mol_ct
for i in range(1,len(pair1)):
#    print i, mol_idx[i],pair1[i],pair2[i],mol_ct
#    raw_input('A:\n\n')
    if (mol_idx[pair1[i]] > -1):
        mol_idx[pair2[i]]= mol_idx[pair1[i]]
    elif (mol_idx[pair2[i]] > -1):
        mol_idx[pair1[i]]= mol_idx[pair2[i]]
    else:
        mol_ct += 1
        mol_idx[pair1[i]]=mol_ct
        mol_idx[pair2[i]]=mol_ct
# Generate list of the atoms for each molecule        
mol_seq = [[]]
for i in range(1,mol_ct):
    mol_seq.append([])
j = mol_idx[pair1[i]]-1
for i in range(l2_num):
    j = mol_idx[i]-1
    mol_seq[j].append(i)
#
# To reconstruct calculate the pair distance if they are on 
# the same molecule and, if needed, add lattice vector
#    
# Need to ident five and six membered rings
#print bond1,bond2  # springs constant and equilbrium lengths
#print len(bdtypes), bdtypes
#
# Now to take the subset of force constants and match them directly
# Lennard-Jones is to be first Note that we can't use the exact hard core as given in hoomd
# We also need to move between Angstroms and kcal/mol to nm and kJ/mol
# Moreover Allinger's MM3 appears to be in terms of r_min and not sigma using the exp-6 Buckingham form 
# sigma = r_min*0.8909 or 1/2**(1./6.)
lj_pair_1 =[]
lj_pair_2 =[]
epsilon = []
sigma = []
dist = 3.0
#k = 0
sig_scale=np.power(2.,-1./6.)  # convert r_min to sigma
for i in range(len(vdw1)):
    for j in range(len(vdw1)):
        lj_pair_1.append("%s" % ptypes[i])
        lj_pair_2.append("%s" % ptypes[j])
#        lj_pair_1.append(aatypes[i])
#        lj_pair_2.append(aatypes[j])
        epsilon_tmp = np.sqrt(vdw2[i]*vdw2[j])   # In kcal/mol 
        sigma_tmp = vdw1[i]+vdw1[j] # in Angstroms
        dist = sigma_tmp
#        E_vdw=epsilon_tmp*(184000.*np.exp(-12.0*dist/sigma_tmp)-2.25*np.power((sigma_tmp/dist),6))
        sigma.append(sig_scale*sigma_tmp)
        epsilon.append(1.1194*epsilon_tmp)   # 1.99*2.25/4.0  
#        E_vdw2=4.*epsilon[k]*(np.power(sigma[k]/dist,12)-np.power((sigma[k]/dist),6))
#        k += 1
#        print aatypes[i],aatypes[j],sigma_tmp,dist, epsilon_tmp,E_vdw,E_vdw2
# Note that these units are for Angstroms and kcal/mol NOT nanometers and kJ/mol
#
atm5ring =[]
atm5ringT =[]
atm5idx =[]
for ring in mol.GetSSSR():
    if ((int(ring.Size()) == 5) and ring.IsAromatic() ):
        temp = ring._path
        temp2 = []
        for i in temp:
            atm5ringT.append(i-1)
            temp2.append(i-1)
        atm5ring.append(temp2)
#
bdtypeid_final = []
bdidx1 = []
bdidx1a = []
bdidx2 = []
bdidx2a = []
bd_E=0.0
for i in range(len(bdtypes_tot)):
    if ('bond5_' in bdtypes_tot[i]):
        bdidx1.append(bdtypes_tot[i])
        bdidx1a.append(i)
    elif ('bond_' in bdtypes_tot[i]):
        bdidx2.append(bdtypes_tot[i])
        bdidx2a.append(i)
for i in range(len(bdtypeid)): 
    if ((pair1[i] in atm5ringT) and (pair2[i] in atm5ringT)):
        for j in range(len(bdidx1)):
            if (bdseq[i] in bdidx1[j]):
                break
        bdtypeid_final.append(bdidx1a[j])
    else:
        for j in range(len(bdidx2)):
            if (bdseq[i] in bdidx2[j]):
                break
        bdtypeid_final.append(bdidx2a[j])
    if (i+1 != len(bdtypeid_final)):
        print 'Something is amiss',(i+1),len(bdtypeid_final)
    p1 =  aax[pair1[i]]-aax[pair2[i]]
    p2 =  aay[pair1[i]]-aay[pair2[i]]
    p3 =  aaz[pair1[i]]-aaz[pair2[i]]
    temp = np.sqrt(p1*p1+p2*p2+p3*p3) 
    temp2 = temp-bond2[bdtypeid_final[i]]
    temp3 = 71.94*bond1[bdtypeid_final[i]]*temp2*temp2
#    print i, pair1[i]+1,pair2[i]+1,bond2[bdtypeid_final[i]],temp,temp3
#    raw_input()
    bd_E += temp3
print "Bond Energy:", bd_E, " in kcal/mol"
#print atm5ring
angtypeid_final = []
angidx1 = []
angidx1a = []
angidx2 = []
angidx2a = []
ang_E= 0.
ang_k_scale=71.94  # 0.021914*57.296*57.296
for i in range(len(angtypes_tot)): # Index all of the angle configurations
    if ('angle5_' in angtypes_tot[i]):
        angidx1.append(angtypes_tot[i])
        angidx1a.append(i)
    elif ('angle_' in angtypes_tot[i]):
        angidx2.append(angtypes_tot[i])
        angidx2a.append(i)
for i in range(len(bangle1)): # Assign angles to configuration type
    test1 = ((bangle1[i] in atm5ringT) and (bangle2[i] in atm5ringT) and (bangle3[i] in atm5ringT))
    test2 = False
    for j in range(len(atm5ring)):
        if ((bangle1[i] in atm5ring[j]) and (bangle2[i] in atm5ring[j]) and (bangle3[i] in atm5ring[j]) and test2 == False):
            test2 = True
    if (test1 and test2):
        for j in range(len(angidx1)):
            if (angseq[i] in angidx1[j]):
                break
        angtypeid_final.append(angidx1a[j])
    else:
        for j in range(len(angidx2)):
            if (angseq[i] in angidx2[j]):
                break
        angtypeid_final.append(angidx2a[j])
    if (i+1 != len(angtypeid_final)):
        print 'Something is amiss',(i+1),len(angtypeid_final)
    p1 = [aax[bangle1[i]]-aax[bangle2[i]],aay[bangle1[i]]-aay[bangle2[i]],aaz[bangle1[i]]-aaz[bangle2[i]]]
    p2 = [aax[bangle3[i]]-aax[bangle2[i]],aay[bangle3[i]]-aay[bangle2[i]],aaz[bangle3[i]]-aaz[bangle2[i]]]
    temp = angle_between(p1, p2)
    temp2 = temp - angt0[angtypeid_final[i]]
#    print i, bangle1[i]+1,aax[bangle1[i]],aay[bangle1[i]],aaz[bangle1[i]]
#    print i, bangle2[i]+1,aax[bangle2[i]],aay[bangle2[i]],aaz[bangle2[i]]
#    print i, bangle3[i]+1,aax[bangle3[i]],aay[bangle3[i]],aaz[bangle3[i]]
#    print bangle1[i]+1,bangle2[i]+1,bangle3[i]+1
#    print i, angt0[angtypeid_final[i]]*57.296,temp*57.296,temp2*57.296,0.021914*57.296*57.296*angk[angtypeid_final[i]]*temp2*temp2
#    raw_input()
    ang_E += 71.94*angk[angtypeid_final[i]]*temp2*temp2
print "Angle Energy:", ang_E, " (kcal/mol)"
#
tortypeid_final = []
for i in range(len(tangle1)): # This not efficient
    if (all(n in atm5ring for n in [tangle1[i],tangle2[i],tangle3[i],tangle4[i]])):
        for j in range(len(tortypes_tot)):
            temp = tortypes_tot[j]
            if (('torsion5_' in temp) and (torseq[i] in temp)):
                tortypeid_final.append(j)
    else:
         for j in range(len(tortypes_tot)):
            temp = tortypes_tot[j]
            if (('torsion_' in temp) and (torseq[i] in temp)):
                tortypeid_final.append(j)
    if ((i+1) != len(tortypeid_final)):  # spans across two rings 
         for j in range(len(tortypes_tot)):
            temp = tortypes_tot[j]
            if (('torsion_' in temp) and (torseq[i] in temp)):
                tortypeid_final.append(j)
    if ((i+1) != len(tortypeid_final)): 
        print 'Something is amiss',i,(i+1),len(tortypeid_final)
#
# MM3 does not have improper torsions but uses out of plane bend
# I could not find the MM4 values so I use those given by  Miller, Jones and Jakowski
#
imptypeid = []; imp_atms = []; imptypes = ['CSCCorSCCS','CCCC']
imp_E = 0.; implist = []
imp_k = [3.59,1.22]
implist.append([imptypes[0],imp_k[0],0.000000001])
implist.append([imptypes[1],imp_k[1],0.000000001])
# The assumuption is that these represent pi-conjugated combinationations 
for i in range(len(tangle1)): # Ad hoc for improper torsions
    temp = atype[tangle1[i]]+atype[tangle2[i]]+atype[tangle3[i]]+atype[tangle4[i]]
    if (('H' in temp) == False):
        p1 = [aax[tangle1[i]],aay[tangle1[i]],aaz[tangle1[i]]]
        p2 = [aax[tangle2[i]],aay[tangle2[i]],aaz[tangle2[i]]]
        p3 = [aax[tangle3[i]],aay[tangle3[i]],aaz[tangle3[i]]]
        p4 = [aax[tangle4[i]],aay[tangle4[i]],aaz[tangle4[i]]]
        temp2 = improper_calc_ang(p1, p2, p3, p4)
        if (temp2 >= 0.0 and temp2 < 1.54): 
#            imp_atms.append([tangle1[i],tangle3[i],tangle2[i],tangle4[i]]) # Doesn't work with HOOMD
            imp_atms.append([tangle2[i],tangle1[i],tangle3[i],tangle4[i]])
        else:
            temp2 = improper_calc_ang(p2, p1, p3, p4)
            if (temp2 >= 0.0 and temp2 < 1.54):
                imp_atms.append([tangle1[i],tangle2[i],tangle3[i],tangle4[i]])
            else:
                print "Improper angle problem; try changing indices", i, temp2
                print "This suggestion has only been tested assuming planarity 1-2-3-4 and 1-3-2-4 or 2-1-3-4"
                raw_input()
        if ('S' in temp):
            imptypeid.append(0)
            imp_E += 0.50*temp2*temp2*imp_k[0]
        else:
            imptypeid.append(1)
            imp_E += 0.50*temp2*temp2*imp_k[1]
        j=len(imptypeid)
#        print j,temp2
#        print j,temp,imptypeid[j-1],imp_atms[j-1]
print "Improper Energy:", imp_E, " (kcal/mol)"
#
# Now to create an ad hoc procedure for defining permanent dipoles in dipole-dipole interactions
#
for i in range (len(dipole_types)):
    [d0,d1,d2]=dipole_types[i]
    if ((d1 == 42 and d2 == 2) or (d1 == 2 and d2 == 42)):  # A C-S dipole
        n = bdtypes.index("S42C2")
        CS_sites = []
        for j in range (len(pair1)): # Find all C-S pairs with the same sulfur atom
            if (bdtypeid[j]==n):
                if (abt[pair1[j]] == 42):
                    CS_sites.append([pair1[j],pair2[j]])  # Sulfur first
                else:
                    CS_sites.append([pair2[j],pair1[j]])
    if ((d1 == 5 and d2 == 2) or (d1 == 2 and d2 == 5)):  # A C-S dipole
        n = bdtypes.index("C2H5")
        CH_sites = []
        for j in range (len(pair1)): # Find all C-S pairs with the same sulfur atom
            if (bdtypeid[j]==n):
                if (abt[pair1[j]] == 2):
                    CH_sites.append([pair1[j],pair2[j]])  # Sulfur first
                else:
                    CH_sites.append([pair2[j],pair1[j]])
CS_sites=sorted(CS_sites,key=lambda tup: tup[0])  # assume thiophene; will be combined
CH_sites=sorted(CH_sites,key=lambda tup: tup[0])  # not necessary
S_site, C_site = zip(*CS_sites)
orient = []
for i in range (l2_num):
    orient.append([[0.0],[1.0],[0.0],[0.0]])
j=0
distSC = 0.0
for i in range (len(CS_sites)/2):
    x1 = aax[S_site[j]]
    y1 = aay[S_site[j]]
    z1 = aaz[S_site[j]]
    x2 = 0.5*(aax[C_site[j]]+aax[C_site[j+1]]) - x1
    y2 = 0.5*(aay[C_site[j]]+aay[C_site[j+1]]) - y1
    z2 = 0.5*(aaz[C_site[j]]+aaz[C_site[j+1]]) - z1
    orient[S_site[j]]=[[0.000],[x2],[y2],[z2]] # points to C which is positive
    distSC += np.sqrt(x2*x2+y2*y2+z2*z2)
#    print j,S_site[j],x1,y1,z1,C_site[j],x2,y2,z2,dist
    j += 2
distSC /= float(len(CS_sites))/2.  # Avg. S C distange |p|=q*dist
# For S C the values are  0.9 and 0.5 which I presume is the charge and the fraction position (halfway)
# For simplicity we put a point dipole on the S, negative to C positive 
C_site, H_site = zip(*CH_sites)
distCH = 0.0
for i in range (len(CH_sites)):
    x1 = aax[H_site[i]]
    y1 = aay[H_site[i]]
    z1 = aaz[H_site[i]]
    x2 = aax[C_site[i]]-x1
    y2 = aay[C_site[i]]-y1
    z2 = aaz[C_site[i]]-z1
    orient[H_site[i]]=[[0.000],[x2],[y2],[z2]] # points along bisector and is negative
    distCH += np.sqrt(x2*x2+y2*y2+z2*z2)
distCH /= float(len(CH_sites))  # Avg. S C distange |p|=q*dist
#
dipole_par2=[]       
for i in range (len(dipole_types)):
    [d0,d1,d2]=dipole_types[i]
    if  ((d1 == 42 and d2 == 2) or (d1 == 2 and d2 == 42)): # Assume thiophene sulfur
        [t1,t2]=dipole_par[i]
        dipole_par2.append([1.7*t1*distSC,t2])
    elif ((d1 == 5 and d2 == 2) or (d1 == 2 and d2 == 5)): # Assume all CH pairs
        [t1,t2]=dipole_par[i]
        dipole_par2.append([t1*distCH,t2])
#    dist=np.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))
#    print j,C_site[i],x1,y1,z1,H_site[i],x2,y2,z2,dist
# Encode quaterion and dipole.
# Although a bit cheesy the easiest thing it so put a point dipole on the sulfur and hydrogen atoms.
# Quantitavely there will be a systematic error but adding it to the C source is cumbersome
#
# It may misrepresent the bond between thiophene rings or thiophene and napthylenes
#    
# Now to establish place atoms within a box and then reconstruct individual molecules for display
#
# Given the unit cell, a, b, c, 90, beta, 90 construct coordinates
# This needs to be modified if triclinic
# Convention for monoclinic construction has a-vector parallel to x-axis, 
# and b-vector parallel to the y-axis  
# But for thin films b and c are in the basal plane and thus in the "y-z" plane
# One needs to be consistent and rotate the atoms accordinging
# Here we will use the crystallographic convention to test NaT2 molecule head
# to head against Tinker
# Building on a substrate will be tricky if still monoclinic
#
# Tidy up the location because Tinker's default is a bit ugly
# From the old structure to the new
xyz = np.array([aax,aay,aaz],dtype=float)
xyz = np.transpose(xyz)
#   Recenter chains outside unit cell periphery
toler = np.array([0.03,0.08,0.05])
xyz = olig_recenter(mol_seq,uc,xyz,toler)
#write_txyz_2('wtest.xyz_2',atype,astring,xyz,abt,uc,header)    
#
a =[uc[0],0.,0.]
b =[0.,uc[1],0.]
c =[(uc[2]*np.cos(np.radians(uc[4]))),0.,(uc[2]*np.sin(np.radians(uc[4])))]
#xy=0.; yz=0.; xz=0.
xy,yz,xz = tilt_calc(a,b,c)
Lx=0.5*uc[0]; Ly=0.5*uc[1]; Lz=0.5*uc[2]
# Lx=25.; Ly=25.; Lz=25. # Oversized to eliminate VdW across grain boundaries
wrapx=np.zeros(l2_num); wrapy=np.zeros(l2_num); wrapz=np.zeros(l2_num)
xyz_new=np.zeros([l2_num,3]); xyz_wp=np.zeros([l2_num,3],dtype=int)
for i in range(l2_num):
    zmin=-Lz
    zmax= Lz
#    print i, aax[i],aay[i],aaz[i],xmin,xmax,ymin,ymax,zmin,zmax 
    for j in range(3):
        xyz_new[i,0] = xyz[i,0]+wrapx[i]
        xyz_new[i,1] = xyz[i,1]+wrapy[i]
        xyz_new[i,2] = xyz[i,2]+wrapz[i]
        xmin=-Lx+(xz-xy*yz)*xyz_new[i,2]+xy*xyz_new[i,1]
        xmax= Lx+(xz-xy*yz)*xyz_new[i,2]+xy*xyz_new[i,1]
        ymin=-Ly+yz*xyz_new[i,2]
        ymax= Ly+yz*xyz_new[i,2]
        while ((xyz[i,0]+wrapx[i]) < xmin):
            wrapx[i] += a[0]
            xyz_wp[i,0] -= 1
            if ((xyz[i,0]+wrapx[i]) >= xmin):
                break
        while ( (xyz[i,0]+wrapx[i]) > xmax):
            wrapx[i] -= a[0]
            xyz_wp[i,0] += 1
            if ((xyz[i,0]+wrapx[i]) <= xmax):
                break
        while ((xyz[i,1]+wrapy[i])  < ymin):
            wrapy[i] += b[1]
            xyz_wp[i,1] -= 1
            if ((xyz[i,1]+wrapy[i]) >= ymin):
                break
        while ((xyz[i,1]+wrapy[i]) > ymax):
            wrapy[i] -= b[1]
            xyz_wp[i,1] += 1
            if ( (xyz[i,1]+wrapy[i]) <= ymax):
                break
        while ((xyz[i,2]+wrapz[i]) < zmin):
            wrapx[i] += c[0]
            wrapz[i] += c[2]
            xyz_wp[i,2] -= 1
            if ((xyz[i,2]+wrapz[i]) >= zmin):
                break
        while ((xyz[i,2]+wrapz[i]) > zmax):
            wrapx[i] -= c[0]
            wrapz[i] -= c[2]
            xyz_wp[i,0] += 1
            if ((xyz[i,2]+wrapz[i]) <= zmax):
                break
#
write_txyz_2('wtest.xyz_2',atype,astring,xyz_new,abt,uc,header)    
#write_txyz('wtest.xyz_2',l2_num,atype,astring,xnew,ynew,znew,abt,uc,header)
# Unwrap atoms to double check
#abc = np.array([[uc[0],0.,0.],[0.,uc[1],0.],[(uc[2]*np.cos(np.radians(uc[4]))),0.,(uc[2]*np.sin(np.radians(uc[4])))]])
#xyz_tmp = xyz_new
#xyz_tmp += np.dot(xyz_wp,abc)
#write_txyz_2('wtest.xyz_3',atype,astring,xyz_tmp,abt,uc,header)
"""
# Test code
for i in range(l2_num):
    xmin=-Lx+(xz-xy*yz)*znew[i]+xy*ynew[i]
    xmax= Lx+(xz-xy*yz)*znew[i]+xy*ynew[i]
    ymin=-Ly+yz*znew[i]
    ymax= Ly+yz*znew[i]
    zmin=-Lz
    zmax= Lz
    tx=int((xnew[i]-xmin)*(xmax-xnew[i]))
    ty=int((ynew[i]-ymin)*(ymax-ynew[i]))
    tz=int((znew[i]-zmin)*(zmax-znew[i]))
    if (tx < 0 or ty < 0 or tz < 0):
        print i,xnew[i],ynew[i],znew[i],tx,ty,tz
        print i,xmin,xmax,ymin,ymax,zmin,zmax
"""
#
# Unwrapped atoms are now contiguous    
#
io_hoomd_params('hoomd_pars.dat','w',bdtypes_tot,angtypes_tot,tortypes_tot,lj_pair_1,lj_pair_2,bond1,bond2,
                uc,astring,atype,abt,header,epsilon,sigma,angk,angt0,tor1,tor2,tor3,tor4,dipole_par2,dipole_types,implist,mol_seq)
#io_hoomd_params('hoomd_pars.dat','r',bdtypes_tot,angtypes_tot,tortypes_tot,lj_pair_1,lj_pair_2,bond1,bond2,
#                uc,astring,atype,abt,header,epsilon,sigma,angk,angt0,tor1,tor2,tor3,tor4)
# Now to write and simple xyz file
 
#    print i, aax[i],aay[i],aaz[i]
#    raw_input()
# Now to check and write out txyz file
# Then to reconstruct continuity of the molecules
# generate molecules earlier and then use pairs and distances to reconstruct
# Finally test hoomd against Tinker
# Make forward progress...finally

#s.configuration.step = 1
#s.particles.N=l2_num
#s.configuration.box=[3, 3, 3, 0, 0, 0]
#s.particles.types=['A','B','C']
#s.particles.types=['42','2','5']
#print s.particles.types  42 2 5 instead of S C H
#s.particles.typeid=[0,0,1,1]
#s.particles.position=[[0,0,0],[1,1,1], [-1,-1,-1], [1,-1,-1]]
#s.particles.position=[[0 for x in range(3)] for y in range(l2_num)]
#
#ang2nm=0.10
ang2nm=1.0
s = gsd.hoomd.Snapshot()
s.configuration.step = 1
s.configuration.box=[2.*Lx*ang2nm,2.*Ly*ang2nm,2.*Lz*ang2nm, xy, xz, yz]  # lxly lz xy xz yz tilts, Now in nanometers 
# All particles must be inside the box, need to deal with molecules that span box walls
# Parameter image there for stores the wrapping information
# from this and any initial wrapping deltas one can reconstruct a contiguous molecule
##############################################################################
# This needs to be done
# http://gsd.readthedocs.io/en/latest/schema-hoomd.html#chunk-configuration/box has
# mathematical transformations explicitly listed
##############################################################################
s.particles.N = l2_num
s.bonds.N = len(bdtypeid_final)
s.bonds.types = bdtypes_tot
s.bonds.typeid = bdtypeid_final
s.angles.N = len(angtypeid_final)
s.angles.types = angtypes_tot
s.angles.typeid = angtypeid_final
s.dihedrals.N = len(tortypeid)
s.dihedrals.types = tortypes_tot
s.dihedrals.typeid = tortypeid_final
s.impropers.N = len(imptypeid)
s.impropers.types = imptypes
s.impropers.typeid = imptypeid
s.bonds.group = [(pair1[i],pair2[i]) for i in range(s.bonds.N)]
s.angles.group = [(bangle1[i],bangle2[i],bangle3[i]) for i in range(s.angles.N)]
s.dihedrals.group = [(tangle1[i],tangle2[i],tangle3[i],tangle4[i]) for i in range(s.dihedrals.N)]
s.impropers.group = [imp_atms[i] for i in range(s.impropers.N)]
s.particles.typeid = [tid[i] for i in range(l2_num)]
s.particles.mass = [(amass[tid[i]]) for i in range(l2_num)]
s.particles.orientation =[orient[i] for i in range(l2_num)] 
#s.particles.types = ['A','B','C']
s.particles.types = [(str(i)) for i in ptypes]
#name='/home/winokur/hoomd_test/test.gsd'
# Transfer information to hoomd gsd file
# Write out the configuration
s.particles.position = [(aax[i]*ang2nm,aay[i]*ang2nm,aaz[i]*ang2nm) for i in range(l2_num)]
u = gsd.hoomd.open('wtest.gsd', mode='wb') # No wrapping for periodic boundary conditions
u.append(s)
s.particles.position = xyz_new*ang2nm
s.particles.image = xyz_wp
u = gsd.hoomd.open('wtest2.gsd', mode='wb')
u.append(s)
#for i in range(l2_num):
#    print i, aax[i]-axwrap[i],aay[i]-aywrap[i],aaz[i]-azwrap[i]
#write_txyz('wtest.xyz_5',l2_num,atype,astring,axwrap,aywrap,azwrap,abt,uc,header)

