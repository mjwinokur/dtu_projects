#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 1 11:27:47 2017

@author: winokur
"""
import os
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
os.system('cd /user/winokur/hoomd_test/')
from mw_library import read_txyz
from obtest import btd_info
from forcefield_decode_mm import mm3_decode
# Decode a tinker file
mfile = '/home/winokur/hoomd_test/NaT2_2.txyz'
l2_num, aax,aay,aaz,atype, astring, abt, uc = read_txyz(mfile)
if (len(uc) != 0):
    raw_input('Stop and remome txyz file unit cell line, openbabel interpreter breaks otherwise')
# Require tinkey key file for unit cell information
uc = [20.502, 5.989, 8.07, 90., 96.8549, 90.]
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
pair1,pair2,bangle1,bangle2,bangle3,tangle1,tangle2,tangle3,tangle4,mol = btd_info("txyz",mfile)
pair1 = [(pair1[i]-1) for i in range(len(pair1))]
pair2 = [(pair2[i]-1) for i in range(len(pair2))]
bangle1 = [(bangle1[i]-1) for i in range(len(bangle1))] # indexing needs to start at zero
bangle2 = [(bangle2[i]-1) for i in range(len(bangle2))]
bangle3 = [(bangle3[i]-1) for i in range(len(bangle3))]
tangle1 = [(tangle1[i]-1) for i in range(len(tangle1))]
tangle2 = [(tangle2[i]-1) for i in range(len(tangle2))]
tangle3 = [(tangle3[i]-1) for i in range(len(tangle3))]
tangle4 = [(tangle4[i]-1) for i in range(len(tangle4))]
print len(pair1),len(bangle1),len(tangle1)
raw_input()
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
amass,bdtypes_tot,bond1,bond2,angtypes_tot,angk,angt0,tortypes_tot,tor1,tor2,tor3 = mm3_decode(fffile,ptypes,aatypes,bdtypes,angtypes,tortypes)
angt0 = map(np.deg2rad, angt0) # convert degrees to radians
# Need to ident five and six membered rings
#print bond1,bond2  # springs constant and equilbrium lengths
#print len(bdtypes), bdtypes
#
# Now to take the subset of force constants and match them directly
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
#        print ring.Size(), ring.IsAromatic(), ring.GetType(), 
#        if ring.IsAromatic():
#            count += 1


# At the moment there are only three bond types  42 2, 2 2, 2 5
# For the pair list we neet

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
s.configuration.box=[50, 50, 50, 0, 0, 0]  # lxly lz xy xz yz tilts 
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
hoomd.data.boxdim(4.*uc[0],4.*uc[1],4.*uc[2],0.0,0.0,0.0,3)
hoomd.init.read_gsd("/home/winokur/hoomd_test/test.gsd", restart=None, frame=0, time_step=None)
#hoomd.init.read_gsd(filename, restart=None, frame=0, time_step=None)
# All particles must be inside the box, need to deal with molecules that span box walls
# The parameter "image" is there for storing the wrapping information
# from this and any initial wrapping deltas one can reconstruct a contiguous molecule
# This needs to be done
# http://gsd.readthedocs.io/en/latest/schema-hoomd.html#chunk-configuration/box has
# mathematical transformations explicitly listed'
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
