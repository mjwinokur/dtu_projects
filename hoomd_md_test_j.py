#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 1 11:27:47 2017

@author: winokur
"""
import hoomd
from hoomd import md
"""
# Independently Read the file out to double check
t = gsd.hoomd.open(name='test.gsd', mode='rb')
snap=t[0]
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
print snap.configuration.box
print snap.configuration.step
"""
# Test hooomd.md
hoomd.context.initialize()
hoomd.data.boxdim(50.,50.,50., 0., 0., 0.,3)
g = hoomd.init.read_gsd("test.gsd", restart=None, frame=0, time_step=None)
lj_pair_1 = []
lj_pair_2 = []
# Two different attempts
#lj_pair_1.append('A')
#lj_pair_2.append('A')
#print g.particles.types
lj_pair_1.append(g.particles.types)
lj_pair_2.append(g.particles.types)
epsilon = []
epsilon.append(0.1000)
sigma = []
sigma.append(3.5)
#
mypair = hoomd.md.pair
pair_coeff = hoomd.md.pair.coeff()
#print lj_pair_1[0],lj_pair_2[0]
#for i in range(len(epsilon)):
#    pair_coeff.set(lj_pair_1[i],lj_pair_2[i],epsilon=epsilon[i], sigma=sigma[i])   
pair_coeff.set(lj_pair_1[0],lj_pair_2[0],epsilon=epsilon[0], sigma=sigma[0])   
print pair_coeff.values
nl=md.nlist.cell()
lj= mypair.lj(r_cut=10.0,nlist=nl)
gall = hoomd.group.all()
md.integrate.mode_standard(dt=0.0)
md.integrate.nve(group=gall)
hoomd.run(1)


