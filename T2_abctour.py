#import numpy as np
#import math
from __future__ import print_function
import os
import sys  
from mw_library import iofiles
from mw_library import read_txyz
from mw_library import read_txyz_strip_uc 
#from mayavi.mlab import *

def read_log_file(name):
#    print('name: ',name
    with open(name, "r") as file:
        text = file.readlines()

    for line in text:
        line=line.strip()  # get rid of \cr\lf at end of line
#	print(l2_num, line
        mystring = line
#	print(l2_num,' mystring:',mystring 
        newline = ' '.join(mystring.split())
        mylist = newline.split(" ")
#	print(mylist 
        if (mylist[0] == 'Total' and mylist[1] == 'Potential'):
            energy = float(mylist[4])
#        if (Total Potential Energy :
#    raw_input('At read_log_file')
    return energy
# 
def write_key(name,uc,pi_flag):
    f = open(name,'w')
    f.write('# Force Field Selection \n')
    f.write('PARAMETERS        ~/MolecularTools/tinker/params/mm3.prm  \n')
    f.write(' \n')
    f.write('# Crystal Lattice And Periodic Boundary \n')
    f.write('A-AXIS   '+str(uc[0])+'\n') 
    f.write('B-AXIS   '+str(uc[1])+'\n') 
    f.write('C-AXIS   '+str(uc[2])+'\n') 
    f.write('ALPHA   '+str(uc[3])+'\n') 
    f.write('BETA    '+str(uc[4])+'\n') 
    f.write('GAMMA   '+str(uc[5])+'\n') 	
    f.write('\n')
    if (pi_flag == 'T3'): 	
        f.write('PISYSTEM  4 6 8 9 11 13 15 17 18 20 21 22 24 26 27 28 30 32 33 34 36 38 39 40 42 43 45 47 49 51 52 54 56 57 58 59 61 62 63 65 67 70 72 73 75 76 77 79 83 85 88 91 94 95 96 97 100 101 102 103 104 107 108 109')	
    if (pi_flag == 'T2'): 	
        f.write('PISYSTEM  2 3 5 7 8 9 11 12 14 16 18 20 21 23 26 27 29 31 32 33 35 36 38 40 42 44 45 47 50 51 53 55 56 57 59 60 62 64 66 68 69 71 74 75 77 79 80 81 83 84 86 88 90 92 93 95')	
    f.close()
    return
###########################################################
#
# Main
#
###########################################################
#
# Example syntax python ../T2_abctour.py -i NaT3_s_tinker_min_flip.txyz -o energy_abcd.dat -t "b0 5.74 c0 7.79 beta0 98 tourtype abcbeta minimize T"
# Example syntax python ~/dtu_projects/T2_abctour.py -i NaT2mw.txyz -o energy_abcd_pi_bonding.dat -t "b0 6.00 c0 8.10 beta0 98 tourtype abcbeta minimize T pibond T2" > current_T2_pibond.dat
# Example syntax python ../T2_abctour.py -i temp4.txyz -o energy_abcd.dat -t "b0 5.74 c0 7.79 beta0 99 tourtype abcbeta minimize T"
# Example syntax python ../T2_abctour.py -i NaT3_s_tinker_min_flip.txyz -o energy_abg.dat -t "b0 5.74 c0 7.79 beta0 101 tourtype albega"
#
# open base momomer with unit cell, the latter will be the central values
# unless changed by a custom line flag -t unitcell a b c alpha beta gamma
#
txyzname = "NaT2mw.txyz" # Base monomer in tinker format and includes the unit cell properties on line 2
txyzname = "wg_test.txyz" # Base monomer in tinker format and includes the unit cell properties on line 2
keyname = "newfile.key"
home = "$HOME"
# Here it important to write a key file and a xyz (without unit cell info)
#49.8844 20.652 6.014 8.17 98.8549 2401
#
# sweep a,b,c,beta
istep=7
minimize = 'F'
pi_flag = 'none'
#pi_flag = 'T3'
#tour_tp = 'albega'
#ofile = "energy_abg.dat"
tour_tp = 'abcbeta'
ofile = "energy_abcd.dat"
ai = float(int(istep/2))+1.
#input()
astep = 0.050; bstep = 0.025; cstep = 0.025
alpha_step =1.00; beta_step = 1.00; gamma_step = 1.00
a0 = 0.0; b0 = 0.0; c0 = 0.0
alpha0 =0.0;beta0 =0.0;gamma0 =0.0
if (len(sys.argv[1:])== 0):
    print("Using defaults", txyzname)
    print("And with no flips or custom build options")
else: 
    ttxyzname,ofile,text = iofiles(sys.argv[1:])
    if (ttxyzname != ""):
        txyzname = ttxyzname
    print("Input ", txyzname,ofile)
    if (text != ""):
        print("with special:",text)
        alist = text.split()
        j = 0
 #       print('alist:',alist,j
        for i in alist:
            if (alist[j] == 'steps'): # Input key file
                istep = int(alist[j+1])
            if (alist[j] == 'astep'): # An arbitrary space
                astep = float(alist[j+1])
                print('astep: ',astep)
            if (alist[j] == 'bstep'): # An arbitrary space
                bstep = float(alist[j+1])
                print('bstep: ',bstep)
            if (alist[j] == 'cstep'): # An arbitrary space
                tiltc = float(alist[j+1])
                print('cstep: ',cstep)
            if (alist[j] == 'betastep'): # An arbitrary space
                tiltc = float(alist[j+1])
                print('betastep: ',beta_step)
            if (alist[j] == 'a0'): # An arbitrary space
                a0 = float(alist[j+1])
                print('a0: ',a0)
            if (alist[j] == 'b0'): # An arbitrary space
                b0 = float(alist[j+1])
                print('b0: ',b0)
            if (alist[j] == 'c0'): # An arbitrary space
                c0 = float(alist[j+1])
                print('c0: ',c0)
            if (alist[j] == 'beta0'): # An arbitrary space
                beta0 = float(alist[j+1])
                print('beta0: ',beta0)
            if (alist[j] == 'minimize'): # An arbitrary space
                minimize = alist[j+1]
                print('minimize: ',minimize)
            if (alist[j] == 'tourtype'): # An arbitrary space
                tour_tp = alist[j+1]
                print('Tour type: ',tour_tp)
            if (alist[j] == 'pibond'): # An arbitrary space
                pi_flag = alist[j+1]
                print('Pi bonding with: ',pi_flag)
            j += 1
uc = read_txyz_strip_uc(txyzname,'newfile0.xyz')
if (len(uc)== 0):
    raw_input('Stop, this tinker file requires unit cell information')
l2_num, aax,aay,aaz,atype, astring, abt, uc2, header = read_txyz('newfile0.xyz')
if (len(uc2) != 0):
    raw_input('Stop and remove txyz file unit cell line, tinker minimize breaks otherwise')
if (a0 == 0.0):
    a0=uc[0]
if (b0 == 0.0):
    b0=uc[1]
if (c0 == 0.0):
    c0=uc[2]
if (alpha0 == 0.0):
    alpha0=uc[3]
if (beta0 == 0.0):
    beta0=uc[4]
if (gamma0 == 0.0):
    gamma0=uc[5]
amin = a0-ai*astep
bmin = b0-ai*bstep
cmin = c0-ai*cstep
alphamin = alpha0 - ai*alpha_step
betamin =  beta0  - ai*beta_step
gammamin = gamma0 - ai*gamma_step
#aS = []
#bS = []
#cS = []
#betaS = []
energy = []
emin=10000000.
nn=0
if (tour_tp == 'albega'):
    ntotal=(istep)*(istep)*(istep)
    for ia in range(istep):
        uc[3] = alphamin+float(ia)*alpha_step 
        for ib in range(istep):
            uc[4] = betamin+float(ib)*beta_step
            for ic in range (istep):
                uc[5] = gammamin+float(ic)*gamma_step
#		create key file
                write_key('newfile.key',uc,pi_flag)
                command='cp '+'newfile0.xyz'+' newfile.xyz'
#                print('System command:',command)
                os.system(command)
                if (minimize == 'T'):
                    command=home+'/MolecularTools/tinker/bin/minimize newfile.xyz 0.001  > newfileM.log'		
                os.system(command)
                command=home+'/MolecularTools/tinker/bin/analyze newfile E > newfile.log'		
                os.system(command)
                e = read_log_file('newfile.log')
                energy.append(e)
                nn += 1
                if (e < emin):
                    emin=e 
                    line = '###### New Minimum E: '+str(emin)+'kcal/mol ######'     
                    print(line)
                    line = "%9.3f%10.4f%10.4f%10.4f" % (e,uc[3],uc[4],uc[5]) 
                    print('Energy: ',line,' Iter ',nn,' of ',ntotal)
#		raw_input()			
                error = os.system('rm newfile.log*')
                error2 = os.system('rm newfile.xyz*')
#		calculate energy
#               minimize
#               calcualte energy
elif (tour_tp == 'abcbeta'):
    ntotal=(istep*istep)*(istep*istep)
    for ia in range(istep):
        uc[0] = amin+float(ia)*astep 
        for ib in range(istep):
            uc[1] = bmin+float(ib)*bstep
            for ic in range (istep):
                uc[2] = cmin+float(ic)*cstep
                for id in range (istep):
                    uc[4] = betamin+float(id)*beta_step
#		create key file
                    write_key('newfile.key',uc,pi_flag)
                    command='cp '+'newfile0.xyz'+' newfile.xyz'
                    os.system(command)
                    if (minimize == 'T'):
                        command=home+'/MolecularTools/tinker/bin/minimize newfile.xyz 0.001  > newfileM.log'		
#                    print('System command:',command)
                    os.system(command)
                    command=home+'/MolecularTools/tinker/bin/analyze newfile E > newfile.log'		
#                    print('System command:',command)
                    os.system(command)
                    e = read_log_file('newfile.log')
#                    print(e)
#                    raw_input('here')
                    energy.append(e)
                    nn += 1
                    if (e < emin):
                        emin=e 
                        line = '###### New Minimum E: '+str(emin)+'kcal/mol ######'     
                        print(line)
                    line = "%9.3f%10.4f%10.4f%10.4f%10.4f" % (e,uc[0],uc[1],uc[2],uc[4]) 
                    print('Energy: ',line,' Iter ',nn,' of ',ntotal)
#		raw_input()			
                    error = os.system('rm newfile.log*')
                    error2 = os.system('rm newfile.xyz*')
#		calculate energy
#               minimize
#               calcualte energy
print('Now to save the results')
#f = open('NaT2mw_energies_flip_min.dat','w')
ii = 0
f = open(ofile,'w')
if (tour_tp == 'albega'):
    for ia in range(istep):
        uc[3] = alphamin+float(ia)*alpha_step 
        for ib in range(istep):
            uc[4] = betamin+float(ib)*beta_step
            for ic in range (istep):
                uc[5] = gammamin+float(ic)*gamma_step
                mynewstring="%9.3f%10.4f%10.4f%10.4f" %  (uc[3],uc[4],uc[5],energy[ii]) 
                f.write(mynewstring+'\n')
                ii += 1			
if (tour_tp == 'abcbeta'):
    for ia in range(istep):
        uc[0] = amin+float(ia)*astep 
        for ib in range(istep):
            uc[1] = bmin+float(ib)*bstep
            for ic in range (istep):
                uc[2] = cmin+float(ic)*cstep
                for id in range (istep):
                    uc[4] = betamin+float(id)*beta_step
                    mynewstring="%9.3f%10.4f%10.4f%10.4f%10.4f" %  (uc[0],uc[1],uc[2],uc[4],energy[ii]) 
                    f.write(mynewstring+'\n')
                    ii += 1			
f.close()
