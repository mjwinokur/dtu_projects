#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 11:26:31 2017

@author: winokur
"""
import os
#import sys  
#
# A python front end for building a tinker anneal input file and initiating a batch process 
#
# Assumption is that both the xyz and the key files already exist
"""

 Enter Cartesian Coordinate File Name :  
 Enter the Initial and Final Temperatures in Degrees K [1000,0] :  
 Enter the Number of Equilibration Steps [0] :  
 Enter the Number of Cooling Protocol Steps [2000] :  
 Use Linear, Sigmoidal or Exponential Cooling Protocol ([L], S or E) :  
 Enter the Time Step Length in Femtoseconds [1.0] :  
 Enter Time between Dumps in Picoseconds [0.1] :  
 Increase Atomic Weights by a Factor of 10^x [x=0.0] :  

~/dtu_projects/test/NaT2mw.xyz 
300,300 
4000 
2000 
L 
1.0 
0.1 
0.0

../params/mm3.prm  #key file must exist so not necessary
"""
# batch_cmd: qsub ~/dtu_projects/batch_tinker_anneal.py
#
print 'A batch front end for tinker anneal'
#
sim_dir = "~/dtu_projects/test/"
sim_base = "NaT2mw"
sim_xyz = sim_base+".xyz"
sim_log = sim_base+".log"
tinker_exe_dir = "~/MolecularTools/tinker/bin/"
tinker_exe = "anneal"
tinker_cmd = tinker_exe_dir+tinker_exe

f = open('input.dat','w')
f.write(sim_dir+sim_xyz+'\n')
f.write('300,300\n')
f.write('4000\n')
f.write('2000\n')
f.write('L\n')
f.write('1.0\n')
f.write('0.1\n')
f.write('0.0\n')
f.close()
#~/MolecularTools/tinker/bin/anneal < ~/dtu_projects/test/input.dat > NaT2mw.log
tinker_input = sim_dir+'input.dat'
print 'Tinker command'
tinker_full_cmd = tinker_cmd+' < '+tinker_input+' > '+sim_log
print tinker_full_cmd
os.system(tinker_full_cmd)
print 'fini'