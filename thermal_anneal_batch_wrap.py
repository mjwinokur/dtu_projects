#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 15:39:38 2017

@author: winokur
"""

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
def tinker_anneal_syntax(sim_dir,sim_log,tinker_exe_dir,tinker_exe,params):
#
# batch_cmd: qsub ~/dtu_projects/batch_tinker_anneal.py
#
#print 'A batch front end for tinker anneal'
#
#sim_dir = "~/dtu_projects/test/"
#sim_base = "NaT2mw"
#sim_xyz = sim_base+".xyz"
#sim_log = sim_base+".log"
#tinker_exe_dir = "~/MolecularTools/tinker/bin/"
#tinker_exe = "anneal"
    tinker_cmd = tinker_exe_dir+tinker_exe
#    f = open('input.dat','w')
#    f.write(sim_dir+sim_xyz+'\n')
#    f.write('300,300\n')
#    f.write('4000\n')
#    f.write('2000\n')
#    f.write('L\n')
#    f.write('1.0\n')
#    f.write('0.1\n')
#    f.write('0.0\n')
#    f.close()
    tinker_input = sim_dir+'/'+params[0]
    f = open(tinker_input,'w')
    for i in range(1,9):
        f.write(params[i])
    f.close()
#~/MolecularTools/tinker/bin/anneal < ~/dtu_projects/test/input.dat > NaT2mw.log
    tinker_log = sim_dir+'/'+ sim_log
#    print 'Tinker command'
    tinker_full_cmd = tinker_cmd+' < '+tinker_input+' > '+tinker_log
#    print tinker_full_cmd
#os.system(tinker_full_cmd)
    return tinker_full_cmd

import os

# This program is a master Python wrapper that called various Python programs
# build a free standing or bulk NaTx film and the initial a long staged simulation at
# fixed density.  A real question is, of course, how to reconcile the lattice expansion
# Experimentally there is almost no change in the layer spacing.
# However, MM3 free standing films, with static density in the layer change their 
# layer spacing.  More than likely this is an artifact.
# One should keep the layer spacing fixed and adjust the intralayer density.
#
# Steps
#   1. Build an NaTx unflipped lattice, 8x3x3 or 8x6x6 (with both txyz and key files)
#   2. Configure simulation
#   3. launch qbatch submission
hdir=os.getenv("HOME")
def_dir=hdir+'/dtu_projects'
# python ../T2_build.py -i NaT3mw.txyz -o NaT3_833_allflip_gap -t "flip all supercell 8 3 3 gap 25.0 0 0 beta 100"
# No gap
# Need to breath the lattice
itype=4
if (itype==1):
    input_txt =' -i '+def_dir+'/NaT3mw.txyz '
    special = ' -t "flip all supercell 8 3 3 gap 0 0 0 " '
    sub_dir='/thermal_NaT2_AF'
    sim_base ='NaT3_833_AF'
if (itype==2):
    input_txt =' -i '+def_dir+'/NaT3mw.txyz '
    special = ' -t "flip none supercell 8 3 3 gap 0 0 0 " '
    sub_dir='/thermal_NaT3_NF'
    sim_base ='NaT3_833_NF'
if (itype==3):
    input_txt =' -i '+def_dir+'/NaT2mw.txyz '
    special = ' -t "flip all supercell 8 3 3 gap 0 0 0 " '
    sub_dir='/thermal_NaT2_AF_c'
    sim_base ='NaT2_833_AF'
if (itype==4):
    input_txt =' -i '+def_dir+'/NaT2mw.txyz '
    special = ' -t "flip none supercell 8 3 3 gap 0 0 0 " '
    sub_dir='/thermal_NaT2_NF'
    sim_base ='NaT2_833_NF'
sim_xyz = sim_base+".xyz"
sim_log = sim_base+".log"
output_txt =' -o '+def_dir+sub_dir+'/'+sim_base+' '
sim_dir = def_dir+sub_dir
cmd1 = 'mkdir '+sim_dir
cmd2 = 'cd '+sim_dir
full_cmd = "python "+def_dir+'/T2_build.py'+input_txt+output_txt+special
#print(cmd1);print(cmd2);print(full_cmd)
os.system(cmd1)
os.system(cmd2)
os.system(full_cmd+' > ~/thermal_anneal_log')
cmd4 = 'cp '+sim_dir+'/'+sim_base+'.txyz '+sim_dir+'/'+sim_base+'.xyz'
os.system(cmd4)
tinker_exe_dir = hdir+"/MolecularTools/tinker/bin/"
tinker_exe = "anneal"
params = []
params.append('input.dat')
params_base='input.dat'
params.append(sim_dir+'/'+sim_xyz+'\n')
params.append('300,300\n')  #2
params.append('4000\n')
params.append('4000\n')
params.append('L\n')
params.append('1.0\n')
params.append('0.1\n')
params.append('0.0\n')
# Need to loop temperature
temps = ['300','350','400','450','500','550','600','650']
ict = len(temps)
#
cmd5='export PATH=$PATH:$/home/winokur/MolecularTools/ffe:/home/winokur/MolecularTools/tinker/bin'
cmd6='LD_LIBRARY_PATH="/home/winokur/MolecularTools/jre/lib/amd64/server:/opt/intel/lib/intel64:$LD_LIBRARY_PATH"'
cmd7='export LD_LIBRARY_PATH'
cmd=[]
#
for i in range(ict):
    params[0]= params_base+'_'+temps[i]
    params[2]=temps[i]+','+temps[i]+'\n'
#    print('2: ',params[2])
    sim_log = sim_base+".log_"+temps[i]
    tinker_full_cmd = tinker_anneal_syntax(sim_dir,sim_log,tinker_exe_dir,tinker_exe,params)
#    print(tinker_full_cmd)
    cmd.append(tinker_full_cmd)
#    os.system(tinker_full_cmd)
#    raw_input()
#print('fini')
tinker_exe = sim_dir+'/tinker.sh'
f = open(tinker_exe,'w')
f.write('!/bin/sh \n')
f.write(cmd5+'\n')
f.write(cmd6+'\n')
f.write(cmd7+'\n')
for i in range(len(cmd)):
    f.write(cmd[i]+'\n')
f.close()
os.system('chmod 755 '+tinker_exe)
#
print(' Use command:  at -q b now -m  -f tinker.sh') 
#
# at -q b -m now -f /home/winokur/dtu_projects/thermal_anneal_batch_wrap.py > /home/winokur/batch.log
# at -q b -m now -f /home/winokur/script.sh >> /home/winokur/batch.log