#import numpy as np
#import math
import os
import sys  
from mw_library import iofiles
from mw_library import read_txyz
from mw_library import read_txyz_strip_uc 
#from mayavi.mlab import *

def read_log_file(name):
#    print 'name: ',name
    with open(name, "r") as file:
        text = file.readlines()

    for line in text:
	line=line.strip();  # get rid of \cr\lf at end of line
#	print l2_num, line
	mystring = line
#	print l2_num,' mystring:',mystring 
	newline = ' '.join(mystring.split())
	mylist = newline.split(" ")
#	print mylist 
	if (mylist[0] == 'Total' and mylist[1] == 'Potential'):
	    energy = float(mylist[4])
#        if (Total Potential Energy :
#    raw_input('At read_log_file')
    return energy
 

###########################################################
#
# Main
#
###########################################################
# open base momomer with unit cell, the latter will be the central values
# unless changed by a custom line flag -t unitcell a b c alpha beta gamma
#
txyzname = "NaT2mw.txyz" # Base monomer in tinker format and includes the unit cell properties on line 2
txyzname = "wg_test.txyz" # Base monomer in tinker format and includes the unit cell properties on line 2
keyname = "newfile.key"
ofile = "energy_abg.dat"
# Here it important to write a key file and a xyx (without unit cell info)
#49.8844 20.652 6.014 8.17 98.8549 2401
#
# sweep a,b,c,beta
istep=5
minimize = 'F'
ai = float(istep/2)+1
astep = 0.050; bstep = 0.025; cstep = 0.025
alpha_step =1.00; beta_step = 1.00; gamma_step = 1.00
a0 = 0.0; b0 = 0.0; c0 = 0.0
alpha0 =0.0;beta0 =0.0;gamma0 =0.0
if (len(sys.argv[1:])== 0):
    print "Using defaults", txyzname,
    print "And with no flips or custom build options"
else: 
    ttxyzname,ofile,text = iofiles(sys.argv[1:])
    if (ttxyzname != ""):
        txyzname = ttxyzname
    print "Files are", txyzname,ofile
    if (text != ""):
        print "with special:",text
        alist = text.split()
        j = 0
 #       print 'alist:',alist,j
        for i in alist:
            if (alist[j] == 'steps'): # Input key file
                istep = int(alist[j+1])
            if (alist[j] == 'astep'): # An arbitrary space
                astep = float(alist[j+1])
                print 'astep: ',astep
            if (alist[j] == 'bstep'): # An arbitrary space
                bstep = float(alist[j+1])
                print 'bstep: ',bstep
            if (alist[j] == 'cstep'): # An arbitrary space
                tiltc = float(alist[j+1])
                print 'cstep: ',cstep
            if (alist[j] == 'betastep'): # An arbitrary space
                tiltc = float(alist[j+1])
                print 'betastep: ',beta_step
            if (alist[j] == 'a0'): # An arbitrary space
                a0 = float(alist[j+1])
                print 'a0: ',a0
            if (alist[j] == 'b0'): # An arbitrary space
                b0 = float(alist[j+1])
                print 'b0: ',b0
            if (alist[j] == 'c0'): # An arbitrary space
                c0 = float(alist[j+1])
                print 'c0: ',c0
            if (alist[j] == 'beta0'): # An arbitrary space
                beta0 = float(alist[j+1])
                print 'beta0: ',c0
            if (alist[j] == 'minimize'): # An arbitrary space
                minimize = float(alist[j+1])
                print 'minimize: ',minimize
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
ntotal=(istep)*(istep)*(istep)
for ia in range(istep):
    at = alphamin+float(ia)*alpha_step 
    for ib in range(istep):
        bt = betamin+float(ib)*beta_step
        for ic in range (istep):
            gt = gammamin+float(ic)*gamma_step
#		create key file
            f = open('newfile.key','w')
            f.write('# Force Field Selection \n')
            f.write('PARAMETERS        /home/winokur/MolecularTools/ffe/../tinker/params/mm3.prm  \n')
            f.write(' \n')
            f.write('# Crystal Lattice And Periodic Boundary \n')
            f.write('A-AXIS   '+str(uc[0])+'\n') 
            f.write('B-AXIS   '+str(uc[1])+'\n') 
            f.write('C-AXIS   '+str(uc[2])+'\n') 
            f.write('ALPHA   '+str(at)+'\n') 
            f.write('BETA    '+str(bt)+'\n') 
            f.write('GAMMA   '+str(gt)+'\n') 		
            f.close()
            command='cp '+'newfile0.xyz'+' newfile.xyz'
            os.system(command)
            if (minimize == 'T'):
                command='/home/winokur/MolecularTools/tinker/bin/minimize newfile.xyz 0.001  > newfileM.log'		
            os.system(command)
            command='/home/winokur/MolecularTools/tinker/bin/analyze newfile E > newfile.log'		
            os.system(command)
            e = read_log_file('newfile.log')
            energy.append(e)
#            print 'Here! ', nn
#            raw_input()
#		aS.append(at)
#		bS.append(bt)
#		cS.append(ct)
#		betaS.append(betat)
            nn += 1
            if (e < emin):
                emin=e 
                line = '###### New Minimum E: '+str(emin)+'kcal/mol ######'     
                print line
            line = "%9.3f%10.4f%10.4f%10.4f" % (e,at,bt,gt) 
            print 'Energy: ',line,' Iter ',nn,' of ',ntotal
#		raw_input()			
            error = os.system('rm newfile.log*')
            error2 = os.system('rm newfile.xyz*')
#		calculate energy
#               minimize
#               calcualte energy
print 'Now to save the results'
f = open(ofile,'w')
#f = open('NaT2mw_energies_flip_min.dat','w')
ii = 0
for ia in range(istep):
    at = alphamin+float(ia)*alpha_step 
    for ib in range(istep):
        bt = betamin+float(ib)*beta_step
        for ic in range (istep):
            gt = gammamin+float(ic)*gamma_step
            mynewstring="%9.3f%10.4f%10.4f%10.4f" %  (at,bt,gt,energy[ii]) 
            f.write(mynewstring+'\n')
            ii += 1			
f.close()
