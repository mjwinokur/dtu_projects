import numpy as np
import math
#import os
#from mayavi.mlab import *

def q_mult(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
    z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
    return w, x, y, z

def normalize(v, tolerance=0.00001):
    mag2 = sum(n * n for n in v)
    if abs(mag2 - 1.0) > tolerance:
        mag = np.sqrt(mag2)
        v = tuple(n / mag for n in v)
    return v
    
def q_conjugate(q):
    w, x, y, z = q
    return (w, -x, -y, -z)

def axisangle_to_q(v, theta):
    v = normalize(v)
    x, y, z = v
    theta /= 2.
    w = np.cos(theta)
    sine = np.sin(theta)
    x = x * sine
    y = y * sine
    z = z * sine
    return w, x, y, z

def q_to_axisangle(q):
    w, v = q[0], q[1:]
    theta = acos(w) * 2.0
    return normalize(v), theta
    
def qv_mult(q1, v1):
    q2 = (0.0,) + v1
    return q_mult(q_mult(q1, q2), q_conjugate(q1))[1:]
        
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
#    raw_input()
    return energy
 
def read_txyz(name):
# open file to get r_x, r_y and r_z coordinates
# open a base file of the monomer 
    print 'name: ',name
    with open(name, "r") as file:
        text = file.readlines()
    global aax
    global aay
    global aaz
    global atype
    global astring
    global header
    aax = []
    aay = []
    aaz = []
    atype = []
    astring = []
    l2_num = 0
    ioff=1

    for line in text:
	line=line.strip();  # get rid of \cr\lf at end of line
	mystring = line
	l2_num += 1
#	print l2_num, mystring 
	if (l2_num == 1):
	    header = line + '\n'
	else:
	    newline = ' '.join(mystring.split())
#	    string.append(newline)
	    mylist = newline.split(" ")
            alen = len(mylist)
            newstring = "%4s" % (mylist[ioff+4])
	    m = alen - 6
#            print(m,alen, mylist )            
	    for l in range(m):  #now to put in bonding configuration
	        n = 6 + l
                temp = int(mylist[n])
#	        mynewstring = mynewstring+' '+str(a_site)
                tstring = "%4s" % (temp) 
	        newstring=newstring+tstring
#            print(newstring)
#
#	    print mylist 
	    atype.append(mylist[ioff])
	    aax.append(float(mylist[ioff+1]))
	    aay.append(float(mylist[ioff+2]))
	    aaz.append(float(mylist[ioff+3]))
	    astring.append(newstring)
	    i=l2_num-2
#	    print aax[i],aay[i],aaz[i],astring[i]
    file.close() 
    l2_num=l2_num-1
    return l2_num
#############################################################

Main

#############################################################
import hoomd
hoomd.context.initialize("")
hoomd.init.create_lattice(unitcell=hoomd.lattice.sc(a=2.0), n=5)
#hoomd.init.create_lattice(unitcell=hoomd.lattice.sc(a=2.0), n=5)
#init.create_random(box=data.boxdim(L=18, xy=0.1, xz=0.2, yz=0.3), N=1000)

#cite.save(file='cite.bib')

# The goal is to read an mm3 list and assign force field parameters from hoomd import md
my_coeffs = md.pair.coeff();
my_force.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
my_force.pair_coeff.set('A', 'B', epsilon=1.0, sigma=2.0)
my_force.pair_coeff.set('B', 'B', epsilon=2.0, sigma=1.0)

name = 'hoomd_trial.py'
    print 'name: ',name
    f = open(name,'w;')
#    global astring
    f.write(header)
    for i in range(l2_num):
        j = i + 1
# Need to break up things into molecules.
# for molecule 1 it will have one neighbor which specify the radius and spring constant for the list
#  	
	espilon_1=
        mynewstring="%4s%3s%12.6f%12.6f%12.6f" % (str(j),atype[i],aax[i],aay[i],aaz[i])
        mynewstring=mynewstring+astring[i]+'\n'
#	print(j,mynewstring)
#	raw_input()
	f.write(mynewstring)
    f.close() 
