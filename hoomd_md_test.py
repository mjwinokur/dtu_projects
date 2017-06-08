#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 1 11:27:47 2017

@author: winokur
"""
#import os
#import openbabel, re  #used in obtest
import numpy as np
import pickle
import hoomd
import hoomd.md
import gsd.hoomd
#import math
#import gsd
from mw_library import tilt_calc 
#from mw_library import write_txyz
from mw_library import write_txyz_2
from mw_library import unwrap_rewrap
from mw_library import olig_recenter
from hoomd_modules import energy_time_plot
#lj_pair_1 = []; lj_pair_2 = []
#uc = [] ; astring = []; atype = []; abt = []; header = []
#epsilon = []; sigma = []; angk = []; angt0 = []
#tor1 = []; tor2 = []; tor3 = []; tor4 = []; bond1 = []; bond2 = []
action = 'r'
name = 'hoomd_pars.dat'
if (action == 'r'):
    f = open(name,'r;')
    bd_types = pickle.load(f)
    ang_types = pickle.load(f)
    tor_types = pickle.load(f)
    lj_pair_1 = pickle.load(f)
    lj_pair_2 = pickle.load(f)
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
    # This tricky because the dipole may vary as a function of ring size but coding this into the pairs definition is not obvious or done as yet
    dipole_par = pickle.load(f)
    dipole_types = pickle.load(f)
    imp_list = pickle.load(f)
    mol_seq = pickle.load(f)
    f.close()
# Read the gsd file to get the necessary parameters
#
# Running the code in nm breaks things, this must be debugged
#
#nm2ang=10. ; ang2nm = 0.10 # Depends on the gsd file format 
nm2ang=1.0; ang2nm = 1.0
t = gsd.hoomd.open(name='wtest2.gsd', mode='rb')
snap=t[0]
print len(t)
l2_num = snap.particles.N
aax = []; aay = []; aaz = []
for i in range(l2_num):
    [x,y,z] = snap.particles.position[i]*nm2ang
    aax.append(x)
    aay.append(y)
    aaz.append(z)
ptypes = []
for i in lj_pair_1:
    if i not in ptypes:
        ptypes.append(i)
# Now to put in the unit cell details, no lower than monoclinic and axes should be conventional to crystallography
a =[uc[0],0.,0.]; b =[0.,uc[1],0.]; c =[(uc[2]*np.cos(np.radians(uc[4]))),0.,(uc[2]*np.sin(np.radians(uc[4])))]
xy,yz,xz = tilt_calc(a,b,c)
#
"""
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
# Unwrap atoms to double check
#xnew=np.zeros(l2_num);ynew=np.zeros(l2_num);znew=np.zeros(l2_num)
abc = np.array([[uc[0],0.,0.],[0.,uc[1],0.],[(uc[2]*np.cos(np.radians(uc[4]))),0.,(uc[2]*np.sin(np.radians(uc[4])))]])
xyz_new = np.array(snap.particles.position)*nm2ang
xyz_wrap =np.array(snap.particles.image,dtype=float) # image has wrap information
xyz_new += np.dot(xyz_wrap,abc)
write_txyz_2('wtest.xyz_3',atype,astring,xyz_new,abt,uc,header)
"""
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
print snap.impropers.group
print snap.impropers.types
print snap.impropers.typeid
print snap.configuration.box
print snap.configuration.step
"""
#  Each hoomd run only allows for a single initialization and so the only
#  way I can figure out how to have a new lattice repeat is to write a new gsd file here
#
# Expand the unit cell and check the calculation
#
abc_mult = [1.01,1.03,1.03]
if (abc_mult[0] != 1.0 or abc_mult[1] != 1.0 or abc_mult[2] != 1.0):  # must be 1.0 or greater
    xyz_new = np.array(snap.particles.position)*nm2ang
    xyz_wrap =np.array(snap.particles.image,dtype=float)
    xyz_new,xyz_wrap,uc = unwrap_rewrap(uc,abc_mult,xyz_wrap,xyz_new)
    snap.configuration.box[0]=uc[0]*ang2nm    
    snap.configuration.box[1]=uc[1]*ang2nm
    snap.configuration.box[2]=uc[2]*ang2nm
#    print " Box, Lx, Ly, Lz:", snap.configuration.box[0:3]
    xyz_name='wwwtest.xyz'
    write_txyz_2(xyz_name,atype,astring,xyz_new,abt,uc,header)
    snap.particles.position = xyz_new*ang2nm
    snap.particles.image = xyz_wrap
    abc = np.array([[uc[0],0.,0.],[0.,uc[1],0.],[(uc[2]*np.cos(np.radians(uc[4]))),0.,(uc[2]*np.sin(np.radians(uc[4])))]])
    xyz_new = np.array(snap.particles.position)*nm2ang
    xyz_wrap =np.array(snap.particles.image,dtype=float)
    xyz_new += np.dot(xyz_wrap,abc) # unwrap the cell
#    write_txyz_2('wwwtest.xyz_2',atype,astring,xyz_new,abt,uc,header)
#   Recenter chains outside unit cell periphery
    toler = np.array([0.03,-0.05,0.05])
    olig_recenter(mol_seq,uc,xyz_new,toler)
    write_txyz_2('wwwtest.xyz_2',atype,astring,xyz_new,abt,uc,header)    
    gsd_input_file="hoomd_temp.gsd"
    u = gsd.hoomd.open(gsd_input_file, mode='wb')
    u.append(snap)
else:
    gsd_input_file="wtest2.gsd"
#print snap.pairs
# Test hooomd.md
#Lx=25.; Ly=25.; Lz=25.
#xy=0.; yz=0.; xz=0.
#hoomd.init.create_lattice(unitcell=hoomd.lattice.sc(a=2.0, type_name='A'), n=10)
#hoomd.data.boxdim(uc[0],uc[1],uc[2], xy, xz, yz,3)
# box with xy,xz,yz must be in the gsd file
hoomd.context.initialize()
g = hoomd.init.read_gsd(gsd_input_file, restart=None, frame=0, time_step=None) # with the wrapped coordinates
#g.particles.types = ptypes
#hoomd.init.read_gsd(filename, restart=None, frame=0, time_step=None)
# All particles must be inside the box, need to deal with molecules that span box walls
# The parameter "image" is there for storing the wrapping information
# from this and any initial wrapping deltas one can reconstruct a contiguous molecule
# http://gsd.readthedocs.io/en/latest/schema-hoomd.html#chunk-configuration/box has
# mathematical transformations explicitly listed'
md=hoomd.md
nl=md.nlist.cell()
nl2=md.nlist.cell()
#mypair = hoomd.md.pair  # Do not use syntax like this
#myljforce = mypair.lj
#nl.reset_exclusions(exclusions = ['1-2', '1-3', '1-4'])
nl.reset_exclusions(exclusions = ['bond', 'angle','dihedral'])
gall = hoomd.group.all()  # all has a reserved function in Python; I don't know why hoomd templates would choose "all" at the example
#
# NOTE:  All distances are in Angstroms and all energies are in kcal/mol
# To properly get the time step we will need to convert energies to Joules/mol and length to nanometers
#pair_coeff = hoomd.md.pair.coeff()
#print pair_coeff.values
#lj= mypair.lj(r_cut=10.0,nlist=nl)
#temp = hoomd.context.current.system_definition.getParticleData().getNameByType(0)
#for i in range(len(epsilon)):
#    pair_coeff.set(lj_pair_1[i],lj_pair_2[i],epsilon=epsilon[i], sigma=sigma[i])   
#pair_coeff.set('A','A',epsilon=epsilon[0], sigma=sigma[0])
#    
kc2kj=4.184
r_ct = 9.0*ang2nm
for i in range (3):
    if (r_ct > 0.4375*uc[i]*ang2nm):
        r_ct = 0.4375*uc[i]*ang2nm
print "r_ct =",r_ct
lj = hoomd.md.pair.lj(r_cut=r_ct, nlist=nl)
for i in range(len(epsilon)):
    lj.pair_coeff.set(lj_pair_1[i],lj_pair_2[i],epsilon=epsilon[i]*kc2kj, sigma=sigma[i]*ang2nm) 
#
# This is a "todo" because hoomd assumes point particles with dipoles
r_ct_dipole = 20.0*ang2nm
for i in range (3):
    if (r_ct_dipole > 0.4375*uc[i]*ang2nm):
        r_ct_dipole = 0.4375*uc[i]*ang2nm
print "r_ct_dipole =",r_ct_dipole
dipole = hoomd.md.pair.dipole(r_cut=r_ct_dipole, nlist=nl2)
for i in range(len(epsilon)):
    dipole.pair_coeff.set(lj_pair_1[i],lj_pair_2[i],mu=0.0, kappa=0.0) 
# Custom assumptions
for i in range (len(dipole_types)):
    [d0,d1,d2]=dipole_types[i]
    if ((d1 == 42 and d2 == 2) or (d1 == 2 and d2 == 42)):
        [t1,t2] = dipole_par[i]
        dipole.pair_coeff.set('42','42',mu=t1*np.sqrt(kc2kj), kappa=0.0)
#        print i, t1
    elif ((d1 == 5 and d2 == 2) or (d1 == 2 and d2 == 5)):
        [t1,t2] = dipole_par[i]
        dipole.pair_coeff.set('5','5', mu=t1*np.sqrt(kc2kj), kappa=0.0)
#        print i, t1
# Only direct terms but cross terms are probably needed as well 
"""
dipole = hoomd.md.pair.dipole(r_cut=r_ct_dipole, nlist=nl2)
for i in range(len(epsilon)):
    dipole.pair_coeff.set(lj_pair_1[i],lj_pair_2[i],mu=0.0, kappa=0.0) 
for i in range(2):
    temp = dipole_types[i]
    temp2=dipole_par[i]
    dipole.pair_coeff.set(str(temp[1]),str(temp[2]),mu=temp2[0], kappa=0.0) 
"""
# Set the bond length force coefficients
#bd_coeff = hoomd.md.bond.coeff()
#for i in range(len(bd_types)):
#    bd_coeff.set(bd_types[i],k=bond1[i], r0=bond2[i])
mm3_bd_scale = 2.*71.94*kc2kj/(ang2nm*ang2nm)  # should be 2.*71.94*kc2kj but why not....I don't know mine is so small
bd_har = hoomd.md.bond.harmonic()
#for i in range(len(bd_types)):
for i in range(len(bd_types)):
    bd_har.bond_coeff.set(bd_types[i],k=(mm3_bd_scale*bond1[i]), r0=bond2[i]*ang2nm)
#
# The MM3 force-field is asymmetric and so there are additional terms
# Tinker gets the correct values
# gnuplot f(x)=71.94*k*(x-x0)**2.*(1.-2.55*(x-x0)+(7./12.)*2.55*(x-x0)**2.)
# It is hard to know what hoomd does exactly but I don't have a one to one agreement
#
#print bd_har.get_energy(gall)   # should be zero
#
# Set the bond angle force coefficients
#ang_coeff = hoomd.md.angle.coeff()
#for i in range(len(ang_types)):
#    ang_coeff.set(ang_types[i],k=angk[i],t0=angt0[i]) # MM3 0.021914*k_theta*(theta-theta0)^2
mm3_ang_scale = 2.*71.9*kc2kj  # 57.3*57.3*0.0219
ang_har=hoomd.md.angle.harmonic()
for i in range(len(ang_types)):
    ang_har.angle_coeff.set(ang_types[i],k=(mm3_ang_scale*angk[i]),t0=angt0[i])
#
# The MM3 force-field is asymmetric and so the expression is
# t(x)=0.021914*k*(x-q0)**2*(1.-0.014*(x-q0)+5.6*10**-5.*(x-q0)**3.-7.0*10**-7*(x-q0)+9.0*10.**-10*(x-q0)**4.)
#
#print ang_har.get_energy(gall)   # should be zero
# Yet to get dihedral dipole-dipole stretch-bend angle-bend
#
opls_dihed = hoomd.md.dihedral.opls()
#opls_dihed_coeff=hoomd.md.dihedral.coeff()
for i in range(len(tor_types)):
    opls_dihed.dihedral_coeff.set(tor_types[i],k1=tor1[i]*kc2kj,k2=tor2[i]*kc2kj,k3=tor3[i]*kc2kj,k4=tor4[i]*kc2kj)
#print opls_dihed.get_energy(gall)   # should be zero
# Set the improper angle 
# MM3 does not have improper torsions but uses out of plane bend
# I could not find the MM4 values so I use those given by  Miller, Jones and Jakowski
# imp_har=hoomd.md.improper.harmonic() # should work but does not grrrrr
imp_har=hoomd.md.improper.harmonic()
for i in range(len(imp_list)):
    [t0,t1,t2]=imp_list[i]
    imp_har.improper_coeff.set(t0, k = t1*kc2kj, chi=t2)
#print imp_har.get_energy(gall)   # should be zero
time_step=0.00
md.integrate.mode_standard(dt=time_step)
nvt = md.integrate.nvt(group=gall, kT=6.*0.80, tau=0.5)
#nvt = md.integrate.nvt(group=gall, kT=0.001, tau=0.5)
hoomd.run(1)
# k_Boltzmann = 0.00831445986144858 kJ/mol/Kelvin
nsteps=500
ncycles=101
tsteps =[]; lj_E =[]; dipl_E = []; bd_E = []; ang_E = []; opls_E = []; imp_E =[]; 
totalE = []
tsteps.append(0); 
lj_E.append(lj.get_energy(gall)/kc2kj)
dipl_E.append(dipole.get_energy(gall)/kc2kj)
bd_E.append(bd_har.get_energy(gall)/kc2kj)
ang_E.append(ang_har.get_energy(gall)/kc2kj)
opls_E.append(opls_dihed.get_energy(gall)/kc2kj)
imp_E.append(imp_har.get_energy(gall)/kc2kj)
print '   lj:',lj_E[0]
print ' dppl:',dipl_E[0]
print ' bond:',bd_E[0]
print '  ang:',ang_E[0]
print ' dihd:',opls_E[0]
print ' impr:',imp_E[0]
totalE.append(lj_E[0]+ bd_E[0] + ang_E[0] + opls_E[0] + dipl_E[0] + imp_E[0])
print "Total:",totalE[0]
time_step=0.005 # If Angstroms then 1.0 is 0.1 picoseconds, if nanometers, then 1.0 is a picosecond
md.integrate.mode_standard(dt=time_step)
#nvt = md.integrate.nvt(group=gall, kT=0.001, tau=0.5)
#md.integrate.nve(group=gall)
image_ct=4
for i in range(1,ncycles):
    hoomd.run(nsteps)
    lj_E.append(lj.get_energy(gall)/kc2kj)
    dipl_E.append(dipole.get_energy(gall)/kc2kj)
    bd_E.append(bd_har.get_energy(gall)/kc2kj)
    ang_E.append(ang_har.get_energy(gall)/kc2kj)
    opls_E.append(opls_dihed.get_energy(gall)/kc2kj)
    imp_E.append(imp_har.get_energy(gall)/kc2kj)
#md.integrate.nve.disable
    print '   lj:',lj_E[i],' dppl:',dipl_E[i],' bond:',bd_E[i],'  ang:',ang_E[i],' dihd:',opls_E[i],' impr:',imp_E[i]
    temp =lj_E[i]+ bd_E[i] + ang_E[i] + opls_E[i] + dipl_E[i] + imp_E[i]
    totalE.append(temp)
    tsteps.append(i*nsteps)
    print "Total Energy:",totalE[i],'after ',tsteps[i],'steps'
    energy_time_plot(tsteps,lj_E,dipl_E,bd_E,ang_E,opls_E,imp_E,totalE)
    if (int(i/5)*5 == i):
        image_ct += 1
        h = g.take_snapshot()
        xyz_name='wwwtest.xyz_'+str(image_ct)
        uc[0]=h.box.Lx; uc[1]=h.box.Ly; uc[2]=h.box.Lz
        xyz_new = np.array(h.particles.position)*nm2ang
        xyz_wrap =np.array(h.particles.image,dtype=float)
        abc = np.array([[uc[0],0.,0.],[0.,uc[1],0.],[(uc[2]*np.cos(np.radians(uc[4]))),0.,(uc[2]*np.sin(np.radians(uc[4])))]])
        xyz_new += np.dot(xyz_wrap,abc)  # unwrap the atoms in the unit cell
        write_txyz_2(xyz_name,atype,astring,xyz_new,abt,uc,header)
        print "Now at t =",float(i)*nsteps*time_step*0.1," picoseconds and writing file ",xyz_name
        raw_input()
        print "Again"
        
raw_input()
#
#nvt.disable()
#md.integrate.nve(group=gall)

h = g.take_snapshot()
uc[0]=h.box.Lx; uc[1]=h.box.Ly; uc[2]=h.box.Lz
xyz_new = np.array(h.particles.position)*nm2ang
xyz_wrap =np.array(h.particles.image,dtype=float) # h.particles.image contains wrap information
abc = np.array([[uc[0],0.,0.],[0.,uc[1],0.],[(uc[2]*np.cos(np.radians(uc[4]))),0.,(uc[2]*np.sin(np.radians(uc[4])))]])
xyz_new += np.dot(xyz_wrap,abc)  # unwrap the atoms in the unit cell
write_txyz_2('wtest.xyz_4',atype,astring,xyz_new,abt,uc,header)

info = hoomd.context.SimulationContext()
print info.on_gpu()
print info.forces

hoomd.context.options.user = 'winokur'
hoomd.context.options.user

system = hoomd.data.system_data
snapshot = system.take_snapshot
print(lj.forces[1])

hoomd.analyze.log("result.log",quantities=['potential_energy'],period=1,overwrite=True)

# Calculate the Lennard Jones energies and pair distances
lj.get_energy(gall)
it=0
for i in range(96):
    tags =np.array([it], dtype=np.int32)
    tags2=np.array([i], dtype=np.int32)
    dist=np.sqrt((aax[i]-aax[it])*(aax[i]-aax[it])+(aay[i]-aay[it])*(aay[i]-aay[it])+(aaz[i]-aaz[it])*(aaz[i]-aaz[it]))
    if (dist < r_ct):
        print it+1,i+1,lj.compute_energy(tags1=tags, tags2=tags2),dist

#

# Calculate pair distances
it=0
i=1
tags =np.array(0, dtype=np.int32)
tags2=np.array(1, dtype=np.int32)
dist=np.sqrt((aax[i]-aax[it])*(aax[i]-aax[it])+(aay[i]-aay[it])*(aay[i]-aay[it])+(aaz[i]-aaz[it])*(aaz[i]-aaz[it]))

t = gsd.hoomd.open(name='test.gsd', mode='rb')
snap=t[0]
for j in range(len(snap.bonds.group)):
    array = snap.bonds.group[j]
    i = int(array[0])
    it = int(array[1])
#    gp = hoomd.group.tags(i,it) # range
    gp = hoomd.group.tag_list(name="a",tags=[i])
    dist=np.sqrt((aax[i]-aax[it])*(aax[i]-aax[it])+(aay[i]-aay[it])*(aay[i]-aay[it])+(aaz[i]-aaz[it])*(aaz[i]-aaz[it]))
    print j,i+1,it+1,bd_har.get_energy(gp),dist
#