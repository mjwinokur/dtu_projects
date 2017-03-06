#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 11:27:47 2017

@author: winokur
"""
#from gsd.hoomd import Snapshot,open,HOOMDTrajectory
# Without these two cryptic lines below the "import gsd" breaks
from gsd.hoomd import Snapshot
z = Snapshot()
import gsd.pygsd  # line below only gets hoomd and not pygsd 
import gsd
import os
#import openbabel, re  #used in obtest
import numpy as np
#import math
os.system('cd /user/winokur/hoomd_test/')
from mw_library import read_txyz
from obtest import btd_info
from forcefield_decode_mm import mm3_decode
# Encode a tinker file
mname = '/home/winokur/hoomd_test/NaT2_2.txyz'
l2_num, aax,aay,aaz,atype, astring, abt, uc = read_txyz(mname)
#name='/home/winokur/hoomd_test/test.gsd'
# Transfer information to hoomd
# Write out the configuration
t = gsd.hoomd.open('/home/winokur/hoomd_test/test.gsd', mode='wb')
s = gsd.hoomd.Snapshot()
#s.configuration.step = 1
#s.particles.N=l2_num
#s.configuration.box=[3, 3, 3, 0, 0, 0]
#s.particles.types=['A','B','C']
#s.particles.types=['42','2','5']
# Now to find what kind of atoms we have
alink = np.empty((l2_num), dtype=int) 
ptypes=[]
tid=[]
j = 0
for i in abt:
    if i not in ptypes:
        ptypes.append(i)
        alink[i] = j
        j += 1        
for i in abt:
    tid.append(alink[i])
#print s.particles.types  42 2 5 instead of S C H
#s.particles.typeid=[0,0,1,1]
#s.particles.position=[[0,0,0],[1,1,1], [-1,-1,-1], [1,-1,-1]]
#s.particles.position=[[0 for x in range(3)] for y in range(l2_num)]
#
# Need unit cell infomation and pair bonding information obabel   
# And force field information
fffile="/opt/tinker/params/mm3.prm"
amass = mm3_decode(fffile,ptypes)
pair1,pair2,bangle1,bangle2,bangle3,tangle1,tangle2,tangle3,tangle4 = btd_info("txyz","/home/winokur/hoomd_test/NaT2_2.txyz")
print len(pair1),len(bangle1),len(tangle1)

s = gsd.hoomd.Snapshot()
#s.particles.types=['42','2','5'] # Currently taken from tinker in reference to MM3 parameters
s.configuration.step = 1
s.configuration.box=[40, 40, 40, 0, 0, 0]  # lxly lz xy xz yz tilts 
# All particles must be inside the box, need to deal with molecules that span box walls
# Parameter image there for stores the wrapping information
# from this and any initial wrapping deltas one can reconstruct a contiguous molecule
# This needs to be done
# http://gsd.readthedocs.io/en/latest/schema-hoomd.html#chunk-configuration/box has
# mathematical transformations explicitly listed
s.particles.N = l2_num
s.bonds.N = len(pair1)
s.angles.N = len(bangle1)
s.dihedrals.N = len(tangle1)
# Something is weird because everything below is Nx2 and not Nx2 Nx3 and Nx4
s.bonds.group = [(pair1[i],pair2[i]) for i in range(s.bonds.N)]
s.angles.group = [(bangle1[i],bangle2[i],bangle3[i]) for i in range(s.angles.N)]
s.dihedrals.group = [(tangle1[i],tangle2[i],tangle3[i],tangle4[i]) for i in range(s.dihedrals.N)]
#
s.particles.position = [(aax[i],aay[i],aaz[i]) for i in range(l2_num)]
s.particles.typeid = [tid[i] for i in range(l2_num)]
s.particles.mass = [(amass[tid[i]]) for i in range(l2_num)]
s.particles.types = [("%s" % i) for i in ptypes]

t.append(s)
#len(t)
#snap = t[0]
#t.append(s)
#t = gsd.hoomd.open('/home/winokur/hoomd_test/test.gsd', mode='wb')

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
print len(t)
print snap.particles.N
print snap.particles.typeid
print snap.particles.position
print snap.particles.types
print snap.configuration.step
print snap.particles.mass
print snap.bonds.group
print snap.angles.group
print snap.dihedrals.group



for s in t[5:-2]:
    print(s.configuration.step,' ')
    
