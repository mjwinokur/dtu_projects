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
import hoomd
from hoomd import md
from gsd.hoomd import Snapshot
z = Snapshot()
import gsd.pygsd  # line below only gets hoomd and not pygsd 
import gsd
from mw_library import read_txyz
from mw_library import read_txyz_strip_uc 
from mw_library import tilt_calc 
from mw_library import write_txyz
from obtest import btd_info
from forcefield_decode_mm import mm3_decode 
# Decode a tinker file
mfile = 'NaT2_2.txyz'
temp_file ='temp.txyz'
uc = read_txyz_strip_uc(mfile,temp_file)
if (len(uc)== 0):
    raw_input('Stop, this tinker file requires unit cell information')
l2_num, aax,aay,aaz,atype, astring, abt, uc2, header = read_txyz(temp_file)
if (len(uc2) != 0):
    raw_input('Stop and remome txyz file unit cell line, openbabel interpreter breaks otherwise')
#
write_txyz('wtest.xyz',l2_num,atype,astring,aax,aay,aaz,abt,uc,header)
#uc = [20.502, 5.989, 8.07, 90., 96.8549, 90.]
# l2_num number of atoms in positions aax,aay,aaz
# atype is the letter element
# astring is the residual bonding information if needed
# abt is currently the MM3 number that identify the bonding configuration
#
# Need unit cell infomation and pair bonding information obabel   
# And force field information
#
fffile ="/opt/tinker/params/mm3.prm"  # MM force constants
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
amass,bdtypes_tot,bond1,bond2,angtypes_tot,angk,angt0,tortypes_tot,tor1,tor2,tor3,vdw1,vdw2,vdwpr1,vdwpr2 = mm3_decode(fffile,ptypes,aatypes,bdtypes,angtypes,tortypes)
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
for ring in mol.GetSSSR():
    if ((int(ring.Size()) == 5) and ring.IsAromatic() ):
        temp = ring._path
        for i in temp:
            j = i-1
            if (j not in atm5ring):
                atm5ring.append(j)
#
bdtypeid_final = []
bdidx1 = []
bdidx1a = []
bdidx2 = []
bdidx2a = []
for i in range(len(bdtypes_tot)):
    if ('bond5_' in bdtypes_tot[i]):
        bdidx1.append(bdtypes_tot[i])
        bdidx1a.append(i)
    elif ('bond_' in bdtypes_tot[i]):
        bdidx2.append(bdtypes_tot[i])
        bdidx2a.append(i)
for i in range(len(bdtypeid)): 
    if ((pair1[i] in atm5ring) and (pair2[i] in atm5ring)):
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
#print atm5ring
angtypeid_final = []
angidx1 = []
angidx1a = []
angidx2 = []
angidx2a = []
for i in range(len(angtypes_tot)):
    if ('angle5_' in angtypes_tot[i]):
        angidx1.append(angtypes_tot[i])
        angidx1a.append(i)
    elif ('angle_' in angtypes_tot[i]):
        angidx2.append(angtypes_tot[i])
        angidx2a.append(i)
for i in range(len(bangle1)): 
    if ((bangle1[i] in atm5ring) and (bangle2[i] in atm5ring) and (bangle3[i] in atm5ring)):
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
# A better way is to assume continuity at the onset and track the wrap count.
#
a =[uc[0],0.,0.]
b =[0.,uc[1],0.]
c =[(uc[2]*np.cos(np.radians(uc[4]))),0.,(uc[2]*np.sin(np.radians(uc[4])))]
xy,yz,xz = tilt_calc(a,b,c)
xy=0.
yz=0.
xz=0.
Lx=0.5*uc[0]
Ly=0.5*uc[1]
Lz=0.5*uc[2]
Lx=25.
Ly=25.
Lz=25.
wrapx=np.zeros(l2_num)
wrapy=np.zeros(l2_num)
wrapz=np.zeros(l2_num)
axwrap=aax
aywrap=aay
azwrap=aaz
wrap_ct=0
for i in range(l2_num):
    xmin=-Lx+(xz-xy*yz)*aaz[i]+xy*aay[i]
    xmax= Lx+(xz-xy*yz)*aaz[i]+xy*aay[i]
    ymin=-Ly+yz*aaz[i]
    ymax= Ly+yz*aaz[i]
    zmin=-Lz
    zmax=+Lz
 #   print i, aax[i],aay[i],aaz[i]
 #   print xmin,xmax,ymin,ymax,zmin,zmax
 #   raw_input()
    while (axwrap[i] < xmin):
        axwrap[i]=axwrap[i]+a[0]
        wrapx[i] +=  1.0
        wrap_ct +=1
        if (axwrap[i] >= xmin):
            break
    while (axwrap[i] > xmax):
        axwrap[i]=axwrap[i]-a[0]
        wrapx[i] += -1.0
        wrap_ct +=1
        if (axwrap[i] <= xmax):
            break
    while (aywrap[i] < ymin):
        aywrap[i]=aywrap[i]+b[1]
        wrapy[i] +=  1.0
        wrap_ct +=1
        if (aay[i] >= ymin):
            break
    while (aywrap[i] > ymax):
        aywrap[i]=aay[i]-b[1]
        wrapy[i] += -1.0
        wrap_ct +=1
        if (aywrap[i] <= ymax):
            break
    while (azwrap[i] < zmin):
        axwrap[i]=axwrap[i]+c[0]
        azwrap[i]=azwrap[i]+c[2]
        wrapz[i] +=  1.0
        wrap_ct +=1
        if (aax[i] >= zmin):
            break
    while (azwrap[i] > zmax):
        axwrap[i]=axwrap[i]-c[0]
        azwrap[i]=azwrap[i]-c[2]
        wrapz[i] += -1.0
        wrap_ct +=1
        if (aaz[i] <= zmax):
            break
#    if (wrap_ct==2):
#        write_txyz('wtest.xyz_2',l2_num,atype,astring,axwrap,aywrap,azwrap,abt,uc,header)
#    print i, axwrap[i],aywrap[i],azwrap[i],wrap_ct
#    raw_input()
#    
write_txyz('wtest.xyz_2',l2_num,atype,astring,axwrap,aywrap,azwrap,abt,uc,header)
#
xunwrap=aax
yunwrap=aay
zunwrap=aaz
for i in range(l2_num):
    xunwrap[i] = axwrap[i] - wrapx[i] * a[0]
    yunwrap[i] = aywrap[i] - wrapx[i] * a[1]
    zunwrap[i] = azwrap[i] - wrapx[i] * a[2]
    xunwrap[i] = axwrap[i] - wrapy[i] * b[0]
    yunwrap[i] = aywrap[i] - wrapy[i] * b[1]
    zunwrap[i] = azwrap[i] - wrapy[i] * b[2]
    xunwrap[i] = axwrap[i] - wrapz[i] * c[0]
    yunwrap[i] = aywrap[i] - wrapz[i] * c[1]
    zunwrap[i] = azwrap[i] - wrapz[i] * c[2]
write_txyz('wtest.xyz_3',l2_num,atype,astring,xunwrap,yunwrap,zunwrap,abt,uc,header)
#
# Unwrapped atoms are now contiguous    
#
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
s = gsd.hoomd.Snapshot()
s.configuration.step = 1
s.configuration.box=[2.*Lx,2.*Ly,2.*Lz, xy, xz, yz]  # lxly lz xy xz yz tilts 
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
s.bonds.group = [(pair1[i],pair2[i]) for i in range(s.bonds.N)]
s.angles.group = [(bangle1[i],bangle2[i],bangle3[i]) for i in range(s.angles.N)]
s.dihedrals.group = [(tangle1[i],tangle2[i],tangle3[i],tangle4[i]) for i in range(s.dihedrals.N)]
#
s.particles.position = [(aax[i],aay[i],aaz[i]) for i in range(l2_num)]
s.particles.typeid = [tid[i] for i in range(l2_num)]
s.particles.mass = [(amass[tid[i]]) for i in range(l2_num)]
s.particles.types = [("%s" % i) for i in ptypes]

#name='/home/winokur/hoomd_test/test.gsd'
# Transfer information to hoomd gsd file
# Write out the configuration
u = gsd.hoomd.open('/home/winokur/hoomd_test/test.gsd', mode='wb')
u.append(s)
""""
# Some example code from the hoomd sample pages
def create_frame(i):
    s = gsd.hoomd.Snapshot()
    s.configuration.step = i+2
    s.particles.N =96
    s.particles.position = np.random.random(size=(96,3))
    return s
t.extend( (create_frame(i) for i in range(2)))
t = gsd.hoomd.open(name='test.gsd', mode='wb')
t.extend( (create_frame(i) for i in range(10)) )
t.append( create_frame(11) )
"""
# Read the file out to double check
t = gsd.hoomd.open(name='/home/winokur/hoomd_test/test.gsd', mode='rb')
snap=t[0]
"""
print len(t)
print snap.particles.N
print snap.particles.typeid
print snap.particles.position
print snap.particles.types
print snap.configuration.step
print snap.particles.mass
print snap.bonds.group
print snap.bonds.types
print snap.bonds.typeid
print snap.angles.group
print snap.angles.types
print snap.angles.typeid
print snap.dihedrals.group
print snap.dihedrals.types
print snap.dihedrals.typeid
"""
# Test hooomd.md
hoomd.context.initialize("")
#hoomd.init.create_lattice(unitcell=hoomd.lattice.sc(a=2.0, type_name='A'), n=10)
# 
hoomd.data.boxdim(2.*Lx,2.*Ly,2.*Lz, xy, xz, yz,3)
hoomd.init.read_gsd("/home/winokur/hoomd_test/test.gsd", restart=None, frame=0, time_step=None)
#hoomd.init.read_gsd(filename, restart=None, frame=0, time_step=None)
# All particles must be inside the box, need to deal with molecules that span box walls
# The parameter "image" is there for storing the wrapping information
# from this and any initial wrapping deltas one can reconstruct a contiguous molecule
# http://gsd.readthedocs.io/en/latest/schema-hoomd.html#chunk-configuration/box has
# mathematical transformations explicitly listed'
pair_coeff = hoomd.md.pair.coeff()
for i in range(len(epsilon)):
    pair_coeff.set(lj_pair_1[i],lj_pair_2[i],epsilon=epsilon[i], sigma=sigma[i])
# Set the bond angle force coefficients
angle_coeff = hoomd.md.angle.coeff()
for i in range(len(angtypes_tot)):
    angle_coeff.set(angtypes_tot[i],k=angk[i],t0=angt0[i])
#
dihed_coeff = hoomd.md.angle.coeff()
for i in range(len(angtypes_tot)):
    dihed_coeff.set(tortypes_tot[i],k1=tor1[i],k2=tor2[i],k3=tor3[i],k4=0.0)
#    
bond_coeff = hoomd.md.bond.coeff()
for i in range(len(bdtypes_tot)):
    bond_coeff.set(bdtypes_tot[i],k=bond1[i], r0=bond2[i])
