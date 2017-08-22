#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 1 11:27:47 2017

@author: winokur
"""
import pickle
import hoomd
import hoomd.md
import gsd.hoomd
# Independently Read the file out to double check
t = gsd.hoomd.open(name='test.gsd', mode='rb')
snap=t[0]
print len(t)
Npt = int(snap.particles.N)
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
bdtypes_tot = []
angtypes_tot = []
tortypes_tot = []
lj_pair_1 = []
lj_pair_2 = []
uc = []
astring = []
atype = []
abt = []
header = []
epsilon = []
sigma = []
angk = []
angt0 = []
tor1 = []
tor2 = []
tor3 = []
tor4 = []
bond1 = []
bond2 = []
action = 'r'
name = 'hoomd_pars.dat'
if (action == 'r'):
    f = open(name,'r;')
    bd_types = pickle.load(f)
    ang_types = pickle.load(f)
    tor_types = pickle.load(f)
    lj_pair_1 = pickle.load(f)
    lj_pair_2 = pickle.load(f)
#    for i in range(len(lj_pair_1)):
#        lj_pair_1[i]=unicode(lj_pair_1[i])
#        lj_pair_2[i]=unicode(lj_pair_2[i])
    bond1 = pickle.load(f)
    bond2 = pickle.load(f)
    uc = pickle.load(f)
    astring = pickle.load(f)
    atype = pickle.load(f)
    abt = pickle.load(f)
    header = pickle.load(f)
    epsilon = pickle.load(f)
    sigma = pickle.load(f)
    angk = pickle.load(f)
    angt0 = pickle.load(f)
    tor1 = pickle.load(f)
    tor2 = pickle.load(f)
    tor3 = pickle.load(f)
    tor4 = pickle.load(f)
    f.close()


g = hoomd.context.initialize("")
"""
snapshot = hoomd.data.make_snapshot(N=snap.particles.N,
                                    box=hoomd.data.boxdim(Lx=50., Ly=50., Lz=50.),
                                    particle_types = snap.particles.types[:],
                                    bond_types = 'bond_C2C2_C2C2')
"""
snapshot = hoomd.data.make_snapshot(N = Npt,
                                    box=hoomd.data.boxdim(Lx=50., Ly=50., Lz=50.),
                                    particle_types =  snap.particles.types,
                                    bond_types = ['bond_C2C2_C2C2','bond_C2H5_H5C2','bond_C2S42_S42C2','bond5_C2C2_C2C2','bond5_C2S42_S42C2'])

snapshot.particles.position[:] = snap.particles.position
snapshot.particles.typeid[:]= snap.particles.typeid
snapshot.bonds.resize(len(snap.bonds.group))
snapshot.bonds.group[:] = snap.bonds.group
#snapshot.replicate(1,20,20)
#import numpy
#snapshot.particles.velocity[:] = numpy.random.normal(0.0,
#  numpy.sqrt(0.8 / 1.0), [snapshot.particles.N, 3])
hoomd.init.read_snapshot(snapshot)
nl = hoomd.md.nlist.cell()
#nl.reset_exclusions(exclusions = ['bond', 'angle','dihedral'])
nl.reset_exclusions(exclusions = [])

lj = hoomd.md.pair.lj(r_cut=10.0, nlist=nl)
for i in range(len(epsilon)):
    lj.pair_coeff.set(lj_pair_1[i],lj_pair_2[i],epsilon=epsilon[i], sigma=sigma[i])   

gall = hoomd.group.all()
hoomd.md.integrate.mode_standard(dt=0.0)
hoomd.md.integrate.nve(group=gall)
hoomd.run(1)

mypair = hoomd.md.pair
lj= mypair.lj(r_cut=10.0,nlist=nl)
#temp = hoomd.context.current.system_definition.getParticleData().getNameByType(0)
pair_coeff = hoomd.md.pair.coeff()
for i in range(len(epsilon)):
    pair_coeff.set('42','2',epsilon=epsilon[i], sigma=sigma[i])   
#    pair_coeff.set(lj_pair_1[i],lj_pair_2[i],epsilon=epsilon[i], sigma=sigma[i])   


gall = hoomd.group.all()

hoomd.md.integrate.mode_standard(dt=0.0)
#hoomd.run(1)  # hoomd complains but then it runs...
hoomd.md.integrate.nve(group=gall)
hoomd.run(1)



dpd = hoomd.md.pair.dpd(r_cut=4.0, nlist=nl, kT=0.8, seed=1)
for i in range (len(snap.bonds.types)):
      dpd.pair_coeff.set('42', '42', A=25.0, gamma = 1.0)
      #dpd.pair_coeff.set('42', '42', A=25.0, gamma = 1.0)
#dpd.pair_coeff.set('42', '2', A=100.0, gamma = 1.0)
#dpd.pair_coeff.set('2', '2', A=25.0, gamma = 1.0)
nl.reset_exclusions(exclusions = [])
harmonic = hoomd.md.bond.harmonic()
harmonic.bond_coeff.set('polymer', k=100.0, r0=0)
hoomd.md.integrate.mode_standard(dt=0.0)
gall = hoomd.group.all();
hoomd.md.integrate.nve(group=gall)
hoomd.run(1)

system = hoomd.data.system_data
snap = system.take_snapshot(g)