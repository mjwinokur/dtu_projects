#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 21:02:38 2017

@author: winokur
"""
import numpy as np
#import math
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
    theta = arccos(w) * 2.0
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
#    global aax
#    global aay
#    global aaz
#    global atype
#    global astring
#    global header
#    global uc
    aax = []
    aay = []
    aaz = []
    atype = []
    abt = []
    astring = []
    uc = []
    l2_num = -1
    ioff=1
    newstring=""

    for line in text:
        line=line.strip();  # get rid of \cr\lf at end of line
        l2_num += 1
        mystring = line
        newline = ' '.join(mystring.split())
        mylist = newline.split(" ")
        alen = len(mylist)
        if (l2_num == 0):
            header = line + '\n'  # Save the first line if needed
        elif (l2_num == 1 and alen == 6):  # Unit cell information is available
            l2_num -= 1
            uc = [(float(mylist[j])) for j in range(6)]
        elif (l2_num > 0):
#    line=line.strip();  # get rid of \cr\lf at end of line
#	mystring = line
#	l2_num += 1
#	print l2_num, mystring 
#	if (l2_num == 1):
#	    header = line + '\n'
#	else:
#	    newline = ' '.join(mystring.split())
#	    string.append(newline)

	    mylist = newline.split(" ")
            alen = len(mylist)
            newstring = "%4s" % (mylist[ioff+5])
	    m = alen - 7
#            print(m,alen, mylist )            
	    for l in range(m):  #now to put in bonding configuration
	        n = 7 + l
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
	    abt.append(int(mylist[ioff+4]))
        astring.append(newstring)
#	    i=l2_num-2
#	    print aax[i],aay[i],aaz[i],astring[i]
    file.close() 
#    l2_num=l2_num-1
    return l2_num, aax,aay,aaz,atype,astring,abt,uc