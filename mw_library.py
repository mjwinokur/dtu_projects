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
#
##############   Quaternion operations ###############
# down to qv_mult
#
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
    theta = np.arccos(w) * 2.0
    return normalize(v), theta
    
def qv_mult(q1, v1):
    q2 = (0.0,) + v1
    return q_mult(q_mult(q1, q2), q_conjugate(q1))[1:]
########################## Find the lowest energy configuration from a Tinker search       
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
######################### Read in a Tinker xyz or txyz file 
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
#    aax = []
#    aay = []
#    aaz = []
    atype = []
    abt = []
    astring = []
    uc = []
    l2_num = -1
    ioff=1
    newstring=""
#
    for line in text:
        mystring=line.strip();  # get rid of \cr\lf at end of line
        l2_num += 1
        newline = ' '.join(mystring.split())
        mylist = newline.split(" ")
        alen = len(mylist)
#        print l2_num,line
#        print l2_num, alen, mylist
        newstring = ""
        if (l2_num == 0):
            header = line  # Save the first line if needed
            lmax = int(mylist[0])
            aax = np.empty((lmax), dtype=object)
            aay = np.empty((lmax), dtype=object)
            aaz = np.empty((lmax), dtype=object)
        elif (l2_num == 1 and alen == 6):  # Unit cell information is available
            l2_num -= 1
            uc = [(float(mylist[j])) for j in range(6)]
        elif (l2_num == lmax+1):
            break
        elif (l2_num > 0):
#            print 'Here alen:',alen
            if (alen > 6):
                newstring = "%4s" % (mylist[ioff+5])
                m = alen - 7
#               print(m,alen, mylist )
                for l in range(m):  
                    n = 7 + l
                    temp = int(mylist[n])
                    tstring = "%5s" % (temp) 
                    newstring=newstring+tstring
            else:
                newstring = -1
#                   print(newstring)
#
#	        print l2_num,mylist 
            i=l2_num-1
            atype.append(mylist[ioff])
            aax[i]=float(mylist[ioff+1])
            aay[i]=float(mylist[ioff+2])
            aaz[i]=float(mylist[ioff+3])
            abt.append(int(mylist[ioff+4]))
            astring.append(newstring)
#            print l2_num,mylist[ioff],float(mylist[ioff+1]),float(mylist[ioff+2]),float(mylist[ioff+3]),int(mylist[ioff+4]),newstring
#            print i, aax[i],aay[i],aaz[i],astring[i]
    file.close() 
    l2_num=len(aax)
    return l2_num,aax,aay,aaz,atype,astring,abt,uc,header
#
######################### Read in an atom xyz file 
#
def read_xyz(name):
# open file to get r_x, r_y and r_z coordinates
# open a base file of the monomer 
    print 'xyz name: ',name
    with open(name, "r") as file:
        text = file.readlines()
    atype = []
    line_num = 0
    aax=[]
    aay=[]
    aaz=[]
    for line in text:
        mystring=line.strip();  # get rid of \cr\lf at end of line
        line_num += 1
        newline = ' '.join(mystring.split())
        mylist = newline.split(" ")
        if (line_num == 1):
            l2_num = int(mylist[0])
        elif (line_num > 2):
            atype.append(mylist[0])
            aax.append(float(mylist[1]))
            aay.append(float(mylist[2]))
            aaz.append(float(mylist[3]))
    file.close() 
#    l2_num=l2_num-1
    return l2_num,aax,aay,aaz,atype
#
def write_xyz(name,atype,aax,aay,aaz):
# open file to get r_x, r_y and r_z coordinates
# open a base file of the monomer 
    print 'xyz write name: ',name
    l2_num = len(atype)
    mynewstring = str(l2_num)+'\n'
    f = open(name,'w;')
    f.write(mynewstring)
    mynewstring = 'cristobalite.xyz'+'\n'
    f.write(mynewstring)
    for i in range(l2_num):
        mynewstring="%6s%10.4f%10.4f%10.4f" % (atype[i],aax[i],aay[i],aaz[i])
        mynewstring=mynewstring+'\n'	
        f.write(mynewstring)
    f.close() 
    return
#
# open a txyz file, extract unit cell parameters, then save a 2nd file without them
#
def read_txyz_strip_uc(name_in,name_out):
    print 'name: ',name_in
    f = open(name_out,'w')
    l2_num = -1
    with open(name_in, "r") as file:
        text = file.readlines()
    uc = []
    for line in text:
#        print line
        mystring = line.strip();  # get rid of \cr\lf at end of line
        l2_num += 1
        newline = ' '.join(mystring.split())
        mylist = newline.split(" ")
        alen = len(mylist)
#        print mylist
        if (l2_num == 0):
            f.write(line)  # write out first line
        elif (l2_num == 1 and alen == 6):  # Unit cell information is available
            l2_num -= 1
            for i in range(6):
                uc.append(float(mylist[i]))
        elif (l2_num > 0):
            f.write(line)
    file.close() 
    f.close()
#    l2_num=l2_num-1
    return uc
#
#
# open a txyz file, extract unit cell parameters, then save a 2nd file without them
#
def read_txyz_info_uc(name_in):
    print 'name: ',name_in
    l2_num = -1
    with open(name_in, "r") as file:
        text = file.readlines()
    uc = []
    for line in text:
#        print line
        mystring = line.strip();  # get rid of \cr\lf at end of line
        l2_num += 1
        newline = ' '.join(mystring.split())
        mylist = newline.split(" ")
        alen = len(mylist)
#        print mylist
        if (l2_num == 1 and alen == 6):  # Unit cell information is available
            l2_num -= 1
            for i in range(6):
                uc.append(float(mylist[i]))
        elif (l2_num > 1):
            break
    file.close() 
#    l2_num=l2_num-1
    return uc
#
def write_tinker_key(name,uc):
    f = open(name,'w')
    f.write('# Force Field Selection \n')
    f.write('PARAMETERS        /home/winokur/MolecularTools/ffe/../tinker/params/mm3.prm  \n')
    f.write(' \n')
    f.write('# Crystal Lattice And Periodic Boundary \n')
    f.write('A-AXIS   '+str(uc[0])+'\n') 
    f.write('B-AXIS   '+str(uc[1])+'\n') 
    f.write('C-AXIS   '+str(uc[2])+'\n') 
    f.write('ALPHA   '+str(uc[3])+'\n') 
    f.write('BETA    '+str(uc[4])+'\n') 
    f.write('GAMMA   '+str(uc[5])+'\n') 		
    f.close()
    return

#
# From http://hoomd-blue.readthedocs.io/en/stable/box.html
# boxMatrix contains an arbitrarily oriented right-handed box matrix.
#
def tilt_calc(a,b,c):
    Lx = np.sqrt(np.dot(a, a))
    a2x = np.dot(a, b) / Lx
    Ly = np.sqrt(np.dot(b,b) - a2x*a2x)
    xy = a2x / Ly
    v0xv1 = np.cross(a, b)
    v0xv1mag = np.sqrt(np.dot(v0xv1, v0xv1))
    Lz = np.dot(c, v0xv1) / v0xv1mag
    a3x = np.dot(a, c) / Lx
    xz = a3x / Lz
    yz = (np.dot(b,c) - a2x*a3x) / (Ly*Lz)
    return xy,yz,xz

def write_txyz(name,l2_num,atype,astring,aax,aay,aaz,abt,uc,header):
# open file to get r_x, r_y and r_z coordinates
# open a base file of the monomer 
    print 'name: ',name
    f = open(name,'w;')
    mynewstring =  str(l2_num)+' '+header
    f.write(mynewstring)
    mynewstring =  "%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f" % (uc[0],uc[1],uc[2],uc[3],uc[4],uc[5])
    f.write(mynewstring+'\n')
    for i in range(l2_num):
        mynewstring="%4s%4s%10.4f%10.4f%10.4f%4s" % (str(i+1),atype[i],aax[i],aay[i],aaz[i],abt[i])
        mynewstring=mynewstring+' '+astring[i]+'\n'
#	print(j,mynewstring)
#	raw_input()
	f.write(mynewstring)
    f.close() 
    return
def write_txyz_2(name,atype,astring,xyz,abt,uc,header):
# open file to get r_x, r_y and r_z coordinates
# open a base file of the monomer 
    print 'name: ',name
    f = open(name,'w;')
    l2_num = len(xyz)
    mynewstring =  str(l2_num)+' '+header
    f.write(mynewstring)
    mynewstring =  "%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f" % (uc[0],uc[1],uc[2],uc[3],uc[4],uc[5])
    f.write(mynewstring+'\n')
    for i in range(l2_num):
        [x,y,z]=xyz[i]
        mynewstring="%4s%4s%10.4f%10.4f%10.4f%4s" % (str(i+1),atype[i],x,y,z,abt[i])
        mynewstring=mynewstring+' '+astring[i]+'\n'
#	print(j,mynewstring)
#	raw_input()
	f.write(mynewstring)
    f.close() 
    return
#
def io_hoomd_params(name,action,bd_types,ang_types,tor_types,lj_pair_1,lj_pair_2,bond1,bond2,uc,astring,atype,abt,header,epsilon,sigma,angk,angt0,tor1,tor2,tor3,tor4,dipole_par,dipole_types,improper,mol_seq):
# open file to get r_x, r_y and r_z coordinates
# open a base file of the monomer 
    import pickle
    print 'name: ',name
    if (action == 'w'):
        f = open(name,'w;')
        pickle.dump(bd_types,f)
        pickle.dump(ang_types,f)
        pickle.dump(tor_types,f)
        pickle.dump(lj_pair_1,f)
        pickle.dump(lj_pair_2,f)
        pickle.dump(bond1,f)
        pickle.dump(bond2,f)
        pickle.dump(uc,f)
        pickle.dump(astring,f)
        pickle.dump(atype,f)
        pickle.dump(abt,f)
        pickle.dump(header,f)
        pickle.dump(epsilon,f)
        pickle.dump(sigma,f)
        pickle.dump(angk,f)
        pickle.dump(angt0,f)
        pickle.dump(tor1,f)
        pickle.dump(tor2,f)
        pickle.dump(tor3,f)
        pickle.dump(tor4,f)
        # This tricky because the dipole may vary as a function of ring size but coding this into the pairs definition is not obvious or done as yet
        pickle.dump(dipole_par,f)
        pickle.dump(dipole_types,f)
        pickle.dump(improper,f)
        pickle.dump(mol_seq,f)
        f.close()
    elif (action == 'r'):
        print ' Does not work'
#        f = open(name,'r;')
#        bd_types = pickle.load(f)
#        f.close()
    return 
"""
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
    return name,action,bd_types,ang_types,tor_types,lj_pair_1,lj_pair_2,bond1,bond2,uc,astring,atype,abt,header,epsilon,sigma,angk,angt0,tor1,tor2,tor3,tor4
"""
# This unwraps the atoms in a PBC cell, modifies the lattice repeat (a,b,c only) and then rewraps the atoms
# If one is using hoomd then the positions need to be updated
def unwrap_rewrap(uc,abc_mult,xyz_wrap,xyz):
    l2_num = len(xyz)
#    a =[uc[0],0.,0.]
#    b =[0.,uc[1],0.]
#    c =[(uc[2]*np.cos(np.radians(uc[4]))),0.,(uc[2]*np.sin(np.radians(uc[4])))]
    abc = np.array([[uc[0],0.,0.],[0.,uc[1],0.],[(uc[2]*np.cos(np.radians(uc[4]))),0.,(uc[2]*np.sin(np.radians(uc[4])))]])
#    xyz_new = np.array(snap.particles.position)*nm2ang
#    xyz_wrap =np.array(snap.particles.image,dtype=float)
    xyz += np.dot(xyz_wrap,abc)
# Unwrapping
#    for k in range(l2_num):
#        aa=float(xwp[k]);bb=float(ywp[k]);cc=float(zwp[k]);
#        aax[k] += aa*a[0]+bb*b[0]+cc*c[0]
#        aay[k] += aa*a[1]+bb*b[1]+cc*c[1]
#        aaz[k] += aa*a[2]+bb*b[2]+cc*c[2]
# 
    uc[0]=uc[0]*abc_mult[0];uc[1]=uc[1]*abc_mult[1];uc[2]=uc[2]*abc_mult[2]
    a =[uc[0],0.,0.]
    b =[0.,uc[1],0.]
    c =[(uc[2]*np.cos(np.radians(uc[4]))),0.,(uc[2]*np.sin(np.radians(uc[4])))]
    xy,yz,xz = tilt_calc(a,b,c)
# Rewrap
    wrapx=np.zeros(l2_num); wrapy=np.zeros(l2_num); wrapz=np.zeros(l2_num)
    xyz_wp=np.zeros([l2_num,3],dtype=int)
    xyz_new=np.zeros([l2_num,3])
    Lx = uc[0]*0.5; Ly = uc[1]*0.5; Lz = uc[2]*0.5;
    for i in range(l2_num):
        zmin=-Lz
        zmax= Lz
#    print i, aax[i],aay[i],aaz[i],xmin,xmax,ymin,ymax,zmin,zmax 
        for j in range(3):
            xyz_new[i,0] = xyz[i,0]+wrapx[i]
            xyz_new[i,1] = xyz[i,1]+wrapy[i]
            xyz_new[i,2] = xyz[i,2]+wrapz[i]
            xmin=-Lx+(xz-xy*yz)*xyz_new[i,2]+xy*xyz_new[i,1]
            xmax= Lx+(xz-xy*yz)*xyz_new[i,2]+xy*xyz_new[i,1]
            ymin=-Ly+yz*xyz_new[i,2]
            ymax= Ly+yz*xyz_new[i,2]
            while ((xyz[i,0]+wrapx[i]) < xmin):
                wrapx[i] += a[0]
                xyz_wp[i,0] -= 1
                if ((xyz[i,0]+wrapx[i]) >= xmin):
                    break
            while ( (xyz[i,0]+wrapx[i]) > xmax):
                wrapx[i] -= a[0]
                xyz_wp[i,0] += 1
                if ((xyz[i,0]+wrapx[i]) <= xmax):
                    break
            while ((xyz[i,1]+wrapy[i])  < ymin):
                wrapy[i] += b[1]
                xyz_wp[i,1] -= 1
                if ((xyz[i,1]+wrapy[i]) >= ymin):
                    break
            while ((xyz[i,1]+wrapy[i]) > ymax):
                wrapy[i] -= b[1]
                xyz_wp[i,1] += 1
                if ((xyz[i,1]+wrapy[i]) <= ymax):
                    break
            while ((xyz[i,2]+wrapz[i]) < zmin):
                wrapx[i] += c[0]
                wrapz[i] += c[2]
                xyz_wp[i,2] -= 1
                if ((xyz[i,2]+wrapz[i]) >= zmin):
                    break
            while ((xyz[i,2]+wrapz[i]) > zmax):
                wrapx[i] -= c[0]
                wrapz[i] -= c[2]
                xyz_wp[i,2] += 1
                if ((xyz[i,2]+wrapz[i]) <= zmax):
                    break
    return xyz_new,xyz_wp,uc
#
def abc_to_rlp(uc):
    import numpy as np
    # Given a, b, c, alpha, beta and gamma we need to calculate a*, b*, c*
#  a . b X c
#  a* = 2 pi/ V_uc (b X c)
#  b* = 2 pi/ V_uc (c X a)
#  c* = 2 pi/ V_uc (a X b)
# G = h a* + k b* + l c*
#  n lambda = 2 d_hkl sin theta
# lamda k / 4 pi = d sin theta 
    cosa=np.cos(np.radians(uc[3])) #; sina=np.sin(np.radians(uc[3]))
    cosb=np.cos(np.radians(uc[4])) #; sinb=np.sin(np.radians(uc[4]))
    cosg=np.cos(np.radians(uc[5])); sing=np.sin(np.radians(uc[5]))
    Vdivab=float(uc[2])*np.sqrt(1.-cosa*cosa-cosb*cosb-cosg*cosg+2.*cosa*cosb*cosg)
#a 0=[uc[0],0.,0.]
#b =[uc[1]*cosg,uc[1]*sing,0.]
#c =[(uc[2]*cosb),uc[2]*(cosa-cosb*cosg)/sing,Vdivab/(sing)]
    ax = uc[0]
    ay = 0.000
    az = 0.000
    bx = uc[1]*cosg
    by = uc[1]*sing
    bz = 0.000
    cx = uc[2]*cosb
    cy = uc[2]*(cosa-cosb*cosg)/sing
    cz = Vdivab/sing
#    cx =-1.*uc[2]*np.sin(np.deg2rad(uc[4]-90.))
#    cy = 0.000
#    cz = uc[2]*np.cos(np.deg2rad(uc[4]-90.))
    a = np.array([ax,ay,az])
    b = np.array([bx,by,bz])
    c = np.array([cx,cy,cz])
    invV=2.0*np.pi/(np.dot(a,np.cross(b,c)))
    astar = invV*(np.cross(b,c))
    bstar = invV*(np.cross(c,a))
    cstar = invV*(np.cross(a,b))
#print astar,bstar,cstar
#Ghkl=hh*astar+kk*bstar+ll*cstar
#    astarhat=astar/np.sqrt(np.dot(astar,astar))
#    bstarhat=bstar/np.sqrt(np.dot(bstar,bstar))
#    cstarhat=cstar/np.sqrt(np.dot(cstar,cstar))
#Gperp=np.dot(Ghkl,astarhat)
#Gmag2=np.dot(Ghkl,Ghkl)
#Gpara=np.sqrt(Gmag2-Gperp*Gperp)
#d=2.*pi/np.sqrt(Gmag2)
#print Ghkl,d,Gperp,Gpara
    return astar,bstar,cstar
#
# open an hklfile and read in the pertinent values
#
def read_hklI_file(name_in):
    print 'name: ',name_in
    hkl = []; I=[]
    with open(name_in, "r") as file:
        text = file.readlines()
    for line in text:
#        print line
        mystring = line.strip();  # get rid of \cr\lf at end of line
        newline = ' '.join(mystring.split())
        mylist = newline.split(" ")
        alen = len(mylist)
#        print mylist
        if (mylist[0] != 'h' and mylist[0] != '+' and alen == 6):
            temp = float(mylist[4])
            if (temp != 0.0000):
                hkl.append([float(mylist[0]),float(mylist[1]),float(mylist[2])])
                I.append(temp*float(mylist[5]))
    file.close()
#    l2_num=l2_num-1
    return hkl,I
#
# If a CIS powder intensity file (i.e., from Mercury) then extract the information
#
def get_hkl_list_from_cisfile(hkl_file,astar,bstar,cstar,Gmax):
    hv = [];  kv = [];  lv = []
    d = []; Ghkl = []; mult = []
    dmin=2.*np.pi/Gmax
    with open(hkl_file, "r") as file:
       text = file.readlines()
       l_num = -1
       arr = np.empty((700), dtype=object) # Define an empty array

    for line in text:
	line=line.strip();  # get rid of \cr\lf at end of line
        l_num += 1
#	print l_num, line
	mystring = line
	mynewlist = ' '.join(mystring.split())
	arr[l_num] = mynewlist
    file.close()
    l_num=l_num-1 # discount 1st line
    ii=-1
    for i in range(l_num+1):
        mystring = arr[i+1].split(" ")
#        print mystring[0],mystring[1],mystring[2],mystring[3]
        dd=float(mystring[3])
	if (dd > dmin):
	    ii += 1
            hh=float(mystring[0])
	    hv.append(mystring[0])
            kk=float(mystring[1])
            kv.append(mystring[1])
            ll=float(mystring[2])
            lv.append(mystring[2])
            d.append(dd)
#    f2[i] = float(mystring[4])*float(mystring[5])  # f^2 times multiplicity
            mult.append(float(mystring[5]))
            Ghkl.append(hh*astar+kk*bstar+ll*cstar)
# Now to sequence the indicies before calculating the structure function
    l_num=ii
    sort_index = np.argsort(d)
# Check the sequence
    for j in range(l_num+1):
        i=sort_index[l_num-j]
#   mynewstring="%4s%4s%3s%3s%3s%13.6f" %  (j,i,hv[i],kv[i],lv[i],float(d[i]))
#    print mynewstring
    return hv,kv,lv,d,Ghkl,mult,sort_index
#
# Generate hkl indicies by brute force
#
def generate_hkl(Gmax,astar,bstar,cstar):
    hv = [];  kv = [];  lv = []
    d = []; Ghkl = []; mult = []
    dmin=2.*np.pi/Gmax
    astarlen=np.sqrt(np.dot(astar,astar))
    bstarlen=np.sqrt(np.dot(bstar,bstar))
    cstarlen=np.sqrt(np.dot(cstar,cstar))
    hmax=int(Gmax/astarlen)+1
    kmax=int(Gmax/bstarlen)+1
    lmax=int(Gmax/cstarlen)+1
    print "h_max,k_max,l_max:",hmax,kmax,lmax
    for ih in range(hmax):
        mval=2.000
        if (ih==0):
            mval=1.0
        for ikk in range(2*kmax+1):
	    ik=ikk-kmax
            for ill in range(2*lmax+1):
                il=ill-lmax
                Gval=ih*astar+ik*bstar+il*cstar
                if (ik==0 and ih==0 and il==0):
                    dval=10000.
                else:
                    dval=2.*np.pi/np.sqrt(np.dot(Gval,Gval))
                if (dval > dmin):
                    hv.append(ih)
                    kv.append(ik)
                    lv.append(il)
                    d.append(dval)
                    Ghkl.append(Gval)
                    mult.append(mval)
# Now to sequence the indicies before calculating the structure function
# based of the d-spacing this groups h,k,l indices together so that the print output is easier to follow
    sort_index = np.argsort(d)
    return hv,kv,lv,d,Ghkl,mult,sort_index
#
# Set up atomic form factors
#
def sequence_atom_ff(atype):
    l2_num=len(atype)
    alink = np.zeros(200, dtype=int) # Define an empty array
    iat = 0
    aff = []
    aform = []
    for i in range(l2_num):
        atom = atype[i]
        if (atom =='H'):
            if (alink[1] == 0): # This tracks which ff have been assign
                iat += 1
                alink[1]=iat
                ffline = '0.489918 20.6593 0.262003 7.74039 0.196767 49.5519 0.049879 2.20159 0.001305'
                ffs = ffline.split(" ")
                for j in range(9):
                    ffs[j]=float(ffs[j])
                aff.append([ffs[0],ffs[1],ffs[2],ffs[3],ffs[4],ffs[5],ffs[6],ffs[7],ffs[8]])
            aform.append(alink[1])  # Points to these parameter for H-type atoms
#aff[1] = '0.489918 20.6593 0.262003 7.74039 0.196767 49.5519 0.049879 2.20159 0.001305
#aff[6] = '2.31 20.8439 31.02 10.2075 1.5886 0.5687 0.865 51.6512 0.2156'
#aff[13]= '6.4202 3.0387 1.9002 0.7426 1.5936 31.5472 1.9646 85.0886 1.1151'
#aff[16]= '6.9053 1.4679 5.2034 22.2151 1.4379 0.2536 1.5863 56.172 0.8669'
        if (atom == 'C'):
            if (alink[6] == 0):
                iat += 1
                alink[6]=iat
                ffline = '2.31 20.8439 1.02 10.2075 1.5886 0.5687 0.865 51.6512 0.2156'
                ffs = ffline.split(" ")
                for j in range(9):
                    ffs[j]=float(ffs[j])
                aff.append([ffs[0],ffs[1],ffs[2],ffs[3],ffs[4],ffs[5],ffs[6],ffs[7],ffs[8]])
            aform.append(alink[6])
        if (atom == 'Al'):
            if (alink[13] == 0):
                iat += 1
                alink[13]=iat
                ffline= '6.4202 3.0387 1.9002 0.7426 1.5936 31.5472 1.9646 85.0886 1.1151'
                ffs = ffline.split(" ")
                for j in range(9):
                    ffs[j]=float(ffs[j])
                aff.append([ffs[0],ffs[1],ffs[2],ffs[3],ffs[4],ffs[5],ffs[6],ffs[7],ffs[8]])
            aform.append(alink[13])
        if (atom == 'S'):
            if (alink[16] == 0):
                iat += 1
                alink[16]=iat
                ffline= '6.9053 1.4679 5.2034 22.2151 1.4379 0.2536 1.5863 56.172 0.8669'
                ffs = ffline.split(" ")
                for j in range(9):
                    ffs[j]=float(ffs[j])
                aff.append([ffs[0],ffs[1],ffs[2],ffs[3],ffs[4],ffs[5],ffs[6],ffs[7],ffs[8]])
            aform.append(alink[16])
    return aform,aff
#
import sys, getopt
#
def iofiles(argv):
   inputfile = ''
   outputfile = ''
   text_ctrl = ''
   try:
      print argv
      opts, args = getopt.getopt(argv,"hi:o:t:",["ifile=","ofile=","text="])
   except getopt.GetoptError:
      print 'test.py -i <inputfile> -o <outputfilebase> -t <"special options">'
      sys.exit(2)
   for opt, arg in opts:
       if opt == '-h':
         print 'test.py -i <inputfile> -o <outputfile> -t <"special options">'
         sys.exit()
       elif opt in ("-i", "--ifile"):
         inputfile = str(arg)
       elif opt in ("-o", "--ofile"):
         outputfile = str(arg)
       elif opt in ("-t", "--text"):
         text_ctrl = str(arg)
#   print 'Input file is "', inputfile
#   print 'Output file is "', outputfile
#   print 'Text is "', text_ctrl
   return inputfile,outputfile,text_ctrl
#
def olig_recenter(mol_seq,uc,xyz_new,toler):
#    toler = np.array([0.03,0.03,0.03])
    n_mol = len(mol_seq)
    a =np.array([uc[0],0.,0.])
    b =np.array([0.,uc[1],0.])
    c =np.array([(uc[2]*np.cos(np.radians(uc[4]))),0.,(uc[2]*np.sin(np.radians(uc[4])))])
    a2=np.dot(a,a)
    b2=np.dot(b,b)
    c2=np.dot(c,c)
    for i in range(n_mol):
        xyz_cm = np.zeros([3])
        for j in mol_seq[i]:
            xyz_cm += xyz_new[j]
        xyz_cm /= len(mol_seq[i])
#        print xyz_cm
        atest = np.dot(a,xyz_cm)/a2 
        btest = np.dot(b,xyz_cm)/b2 
        ctest = np.dot(c,xyz_cm)/c2 
        xm = np.rint(atest+toler[0])
        ym = np.rint(btest+toler[1])
        zm = np.rint(ctest+toler[2])
        if (toler[0] != 0. and xm != 0.0):
            for j in mol_seq[i]:
                xyz_new[j] -= xm*a
        if (toler[1] != 0. and ym != 0.0):
            for j in mol_seq[i]:
                xyz_new[j] -= ym*b
        if (toler[2] != 0. and zm != 0.0):
            for j in mol_seq[i]:
                xyz_new[j] -= zm*c
#        if (xm !=0 or ym != 0.0 or zm !=0)
#            print i, atest,btest,ctest,np.rint(atest+toler[0]),np.rint(btest+toler[1]),np.rint(ctest+toler[2])
    for i in range(n_mol):
        xyz_cm = np.zeros([3])
        for j in mol_seq[i]:
            xyz_cm += xyz_new[j]
        xyz_cm /= len(mol_seq[i])
#        print xyz_cm
        atest = np.dot(a,xyz_cm)/a2 
        btest = np.dot(b,xyz_cm)/b2 
        ctest = np.dot(c,xyz_cm)/c2 
        xm = np.rint(atest+toler[0])
        ym = np.rint(btest+toler[1])
        zm = np.rint(ctest+toler[2])
        if (toler[0] != 0. and xm != 0.0):
            for j in mol_seq[i]:
                xyz_new[j] -= xm*a
        if (toler[1] != 0. and ym != 0.0):
            for j in mol_seq[i]:
                xyz_new[j] -= ym*b
        if (toler[2] != 0. and zm != 0.0):
            for j in mol_seq[i]:
                xyz_new[j] -= zm*c
#        print i, btest,np.rint(btest+toler[1])
    return xyz_new
#
def read_data_file(name):
# open file to get r_x, r_y and r_z coordinates
# open a base file of the T2 monomer 
    print 'name: ',name
    with open(name, "r") as file:
        text = file.readlines()
    j=0
    xd=[]
    yd0=[]
    yd1=[]
    yd2=[]
    istep=5
    step=1./(float(istep))
    ict=0

    for line in text:
	line=line.strip();  # get rid of \cr\lf at end of line
#	print l2_num, line
	mystring = line
#	print l2_num,' mystring:',mystring 
	newline = ' '.join(mystring.split())
	mylist = newline.split(" ")
#	print ict,mylist 
        if (j==0):
	    xtemp=float(mylist[0])
	    y0temp=float(mylist[1])
	    y1temp=float(mylist[2])
	    y2temp=float(mylist[3])
	    j = 1
	else:
	    j += 1
	    xtemp=xtemp+float(mylist[0])
	    y0temp=y0temp+float(mylist[1])
	    y1temp=y1temp+float(mylist[2])
	    y2temp=y2temp+float(mylist[3])
	if (j==istep):
	    j = 0
            xd.append(xtemp*step)
            yd0.append(y0temp*step)
            yd1.append(y1temp*step)
            yd2.append(y2temp*step)
#	    print ict,xd[ict],yd0[ict],yd1[ict],yd2[ict]
	    ict += 1

#    g = plt.figure(2)
#    plt.subplot()
#plt = fig.add_subplot(211)
#    g.show
#    raw_input()
    return ict,xd,yd0,yd1,yd2
#
def merge_txyz(txyzname,subtxyzname,ashift,bshift,cshift,special):
    print txyzname
    l2_num,aax,aay,aaz,atype,astring,abt,uc,header = read_txyz(txyzname)
    print subtxyzname
    s2_num,sx,sy,sz,satype,sastring,sabt,suc,header = read_txyz(subtxyzname)
#    raw_input('here')
    ax =[];ay=[];az=[]
    for i in range(l2_num):
        ax.append(aax[i])
        ay.append(aay[i])
        az.append(aaz[i])
    for i in range(s2_num):
        ax.append(sx[i]+ashift)
        ay.append(sy[i]+bshift)
        az.append(sz[i]+cshift)
        abt.append(sabt[i])
        atype.append(satype[i])
        newline = ' '.join(sastring[i].split())
        mylist = newline.split(" ")
        alen = len(mylist)
        for j in range(alen):
            mylist[j]=int(mylist[j])+ l2_num 
        newstring = "%4s" % (str(mylist[0]))       
#        print i,sastring[i]
        for n in range(1,alen):
            tstring ="%5s" % (str(mylist[n]))
            newstring = newstring + tstring
        astring.append(newstring)
#        print 'new string:',newstring
#        raw_input()
    return l2_num,s2_num,ax,ay,az,atype,astring,abt,uc,header
