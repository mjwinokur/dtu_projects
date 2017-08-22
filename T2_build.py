#
#from mw_library import q_mult
#from mw_library import normalize
#from mw_library import q_conjugate
from mw_library import axisangle_to_q
#from mw_library import q_to_axisangle
from mw_library import qv_mult
from mw_library import read_txyz
from mw_library import write_txyz
from mw_library import iofiles
from mw_library import merge_txyz
from rbuildTx import rebuild_Tx
import numpy as np
import sys
import random
#import scipy.optimize
#import functools
#import math
#import os
#from mayavi.mlab import *
#########################################################################
#
# Main:  This program builds a NaT2 supercell from a base unit cell
# in Tinker format and then, using hoomd_gsd_test, takes this information and 
# that from openbabel to code two files for use in hoomd.  The first files, a 
# a gsd format file
# contains the information for the atoms and connectivity and encodes for
# periodic boundary conditions.  The second, a Python dump, contains the 
# interaction coefficients and other key details.   
#
#########################################################################
def plane(x, y, params):
    a = params[0]
    b = params[1]
    c = params[2]
    z = a*x + b*y + c
    return z

def error(params, points):
    result = 0
    for (x,y,z) in points:
        plane_z = plane(x, y, params)
        diff = abs(plane_z - z)
        result += diff**2
    return result
#
# open base momomer as an input 
#txyzname = "NaT3_s_tinker_min.xyz" # From tinker minimization
#txyzname = "NaT3_s_tinker_min.xyz_2" # From tinker minimization
txyzname = "NaT3mw.txyz" # Base monomer in tinker format and includes the unit cell properties on line 2
#
# Two output files,  base name (to the left must be the same)
#
# 
txyznewname = "wf_test.txyz"  # or just ztest.xyz   ########## Look at me!
tnewkey="wf_test.key"  ######### Look at me!
# Repeat of the unit cell along the principle lattice vectors,  currently 3x3x3
lenx=3; leny=3; lenz=3
rough = 'F'
tilta = 0.0; tiltb = 0.0; tiltc = 0.0
alpha=0.; beta=0.; gamma = 0.0
#lenx=6; leny=2; lenz=2
agap = 0.0; bgap = 0.0; cgap = 0.0 # in Angstroms
iflip = 0; jflip = 1; kflip =1;
iflip = 'none'; # No flips needed
substrate = 'none'
nspecial = 0; nsp_flag = 'T'; nsp=[]
#iflip = 'all'
if (len(sys.argv[1:])== 0):
    print "Using defaults", txyzname,txyznewname,tnewkey
    print "With a cell that is ",lenx,"x",leny,"x",lenz
    print "And with no flips or custom build options"
else: 
    ttxyzname,ttxyznewname,text = iofiles(sys.argv[1:])
    if (ttxyznewname != ""):
        txyznewname = ttxyznewname+'.txyz'
        tnewkey = ttxyznewname+'.key'
    if (ttxyzname !=""):
        txyzname = ttxyzname
    print "Files are", txyzname,txyznewname,tnewkey
    if (text != ""):
        print "with special:",text
        alist = text.split()
        j = 0
 #       print 'alist:',alist,j
        for i in alist:
            if (alist[j]== 'flip' and alist[j+1] == 'all'):
                iflip ='all'
            elif (alist[j]== 'flip' and alist[j+1] == 'random'):
                iflip ='random'
            elif (alist[j]== 'flip' and alist[j+1] == 'layers'):
                iflip ='layers'
                jflip = alist[j+2]
                jflip = jflip.split(',')
                for ii in range (len(jflip)):
                    jflip[ii]= int(jflip[ii])
            elif (alist[j]== 'flip' and alist[j+1] == 'none'):
                iflip ='none'
            elif (alist[j]== 'flip'): # Flip one or half the molecules in the unit cell
                iflip = int(alist[j+1])
                jflip = int(alist[j+2])
                kflip = int(alist[j+3])
                print iflip,jflip,kflip
            if (alist[j] == 'supercell'): # Number of unit cells
                lenx = int(alist[j+1])
                leny = int(alist[j+2])
                lenz = int(alist[j+3])
                print lenx,leny,lenz
            if (alist[j] == 'gap'): # An arbitrary space
                agap = float(alist[j+1])
                bgap = float(alist[j+2])
                cgap = float(alist[j+3])
                print 'gaps: ',agap,bgap,cgap
            if (alist[j] == 'tilta'): # An arbitrary space
                tilta = float(alist[j+1])
                print 'tilta: ',tilta
            if (alist[j] == 'tiltb'): # An arbitrary space
                tiltb = float(alist[j+1])
                print 'tiltb: ',tiltb
            if (alist[j] == 'tiltc'): # An arbitrary space
                tiltc = float(alist[j+1])
                print 'tiltc: ',tiltc
            if (alist[j] == 'alpha'): # An arbitrary space
                alpha = float(alist[j+1])
                print 'alpha: ',alpha
            if (alist[j] == 'beta'): # An arbitrary space
                beta = float(alist[j+1])
                print 'beta: ',beta
            if (alist[j] == 'gamma'): # An arbitrary space
                gamma = float(alist[j+1])
                print 'gamma: ',gamma
            if (alist[j] == 'substrate'): # A substrate is to be added
                substrate = alist[j+1]
                ashift = float(alist[j+2])
                bshift = float(alist[j+3])
                cshift = float(alist[j+4])
                special ='none'
                print 'substrate: ',substrate
            j += 1
#
l2_num,aax,aay,aaz,atype,astring,abt,uc,header = read_txyz(txyzname)
if (len(uc) != 6):
    raw_input("Tinker file must have the unit cell parameters")
if (alpha != 0.):
    uc[3]=alpha
if (beta != 0.):
    uc[4]=beta
if (gamma != 0.):
    uc[5]=gamma
nct=lenx*leny*lenz*l2_num
datamean = []
bfline = []
if (l2_num == 96):
    monomer = 'T2'
    xt = np.append(np.asarray(aax[0:24],dtype=float),np.asarray(aax[48:72],dtype=float))
    yt = np.append(np.asarray(aay[0:24],dtype=float),np.asarray(aay[48:72],dtype=float))
    zt = np.append(np.asarray(aaz[0:24],dtype=float),np.asarray(aaz[48:72],dtype=float))
    data = np.concatenate((xt[:,np.newaxis],yt[:,np.newaxis],zt[:,np.newaxis]),axis=1)
# Calculate the mean of the points, i.e. the 'center' of the cloud
    datamean.append(data.mean(axis=0))
# Do an SVD on the mean-centered data.
    uu, dd, vv = np.linalg.svd(data - datamean[0])
    bfline.append(vv[0])
#
    xt = np.append(np.asarray(aax[24:48],dtype=float),np.asarray(aax[72:96],dtype=float))
    yt = np.append(np.asarray(aay[24:48],dtype=float),np.asarray(aay[72:96],dtype=float))
    zt = np.append(np.asarray(aaz[24:48],dtype=float),np.asarray(aaz[72:96],dtype=float))
    data = np.concatenate((xt[:,np.newaxis],yt[:,np.newaxis],zt[:,np.newaxis]),axis=1)
# Calculate the mean of the points, i.e. the 'center' of the cloud
    datamean.append(data.mean(axis=0))
# Do an SVD on the mean-centered data.
    data = data - datamean[0]
    uu, dd, vv = np.linalg.svd(data)
    bfline.append(vv[0])
# vv[0] is the best-fit line
elif (l2_num == 110):
    monomer = 'T3'
    xt = np.asarray(aax[0:55],dtype=float)
    yt = np.asarray(aay[0:55],dtype=float)
    zt = np.asarray(aaz[0:55],dtype=float)
    data = np.concatenate((xt[:,np.newaxis],yt[:,np.newaxis],zt[:,np.newaxis]),axis=1)
# Calculate the mean of the points, i.e. the 'center' of the cloud
    datamean.append(data.mean(axis=0))
# Do an SVD on the mean-centered data.
    uu, dd, vv = np.linalg.svd(data - datamean[0])
    bfline.append(vv[0])
#
    xt = np.asarray(aax[55:110],dtype=float)
    yt = np.asarray(aay[55:110],dtype=float)
    zt = np.asarray(aaz[55:110],dtype=float)
    data = np.concatenate((xt[:,np.newaxis],yt[:,np.newaxis],zt[:,np.newaxis]),axis=1)
# Calculate the mean of the points, i.e. the 'center' of the cloud
    datamean.append(data.mean(axis=0))
# Do an SVD on the mean-centered data.
    data = data - datamean[0]
    uu, dd, vv = np.linalg.svd(data)
    bfline.append(vv[0])  # This is the best-fit line
#
# Vectors perpendicular to the bfline
# Find best-fit plane
#fun = functools.partial(error, points=data)
#params0 = [0, 0, 0]
#res = scipy.optimize.minimize(fun, params0)
#a = res.x[0]; b = res.x[1]; c = res.x[2]
#point  = np.array([0.0, 0.0, c])
#norm_out = np.cross([1,0,a], [0,1,b])  #Normal to the plane
#d = -point.dot(norm_out)
    
# The oligomer orientation may need to be altered
#
nct = 0
if (monomer == 'T2'):
    itt=48
    imap = np.empty(2*itt,dtype=int)
    for i in range(24):  # the atoms of a monomer are not contiguous
        imap[i]=i
        imap[i+24]=i+itt
        imap[i+48]=i+24
        imap[i+72]=i+itt+24
elif (monomer == 'T3'):
    itt=55 # T3 monomer is continuous
    imap = np.empty(2*itt,dtype=int)
    for i in range(2*itt):  # the atoms of a monomer are contiguous
        imap[i]=i
# Use quaternions to rotate molecules
#  With zero this does nothing but for NaT3 it will matter
#q1 = axisangle_to_q([aax[i1]-aax[i2],aay[i1]-aay[i2],aaz[i1]-aaz[i2]],0.0)
q1 = axisangle_to_q(bfline[0],0.0)
# This will rotate one of the two monomer 180 deg to give head to head
#q1R = axisangle_to_q([aax[i1]-aax[i2],aay[i1]-aay[i2],aaz[i1]-aaz[i2]],np.pi)
q1R = axisangle_to_q(bfline[0],np.pi)
aaxR =[]; aayR =[]; aazR =[]
for i in range(l2_num):  # Generate base for flipped array
    aaxR.append(aax[i]); aayR.append(aay[i]); aazR.append(aaz[i])
#print q1
#v2 = qv_mult(q1,(xt[1],yt[1],zt[1])) 
#print [xt[1],yt[1],zt[1]],v2
####### Needs modification if triclinic unit cells are to be used
cosa=np.cos(np.radians(uc[3])) #; sina=np.sin(np.radians(uc[3]))
cosb=np.cos(np.radians(uc[4])) #; sinb=np.sin(np.radians(uc[4]))
cosg=np.cos(np.radians(uc[5])); sing=np.sin(np.radians(uc[5]))
Vdivab=uc[2]*np.sqrt(1.-cosa*cosa-cosb*cosb-cosg*cosg+2.*cosa*cosb*cosg)
a =[uc[0],0.,0.]
b =[uc[1]*cosg,uc[1]*sing,0.]
c =[(uc[2]*cosb),uc[2]*(cosa-cosb*cosg)/sing,Vdivab/(sing)]
#c =[(uc[2]*np.cos(np.radians(uc[4]))),0.,(uc[2]*np.sin(np.radians(uc[4])))]
[dmx,dmy,dmz] = datamean[0]
#for j in range(l2_num):  
for j in range(itt):  
    i=imap[j]  # This flips just the 2nd monomer in the unit cell
    xt=aax[i]-dmx;yt=aay[i]-dmy;zt=aaz[i]-dmz
#    [xt,yt,zt]=data[i]-datamean[0]
#    v2 = qv_mult(q1,(xt,yt,zt))  #No flip
#    aax[i]=v2[0]; aay[i]=v2[1]; aaz[i]=v2[2]
    v2R = qv_mult(q1R,(xt,yt,zt)) #Flip
    aaxR[i]=v2R[0]+dmx; aayR[i]=v2R[1]+dmy; aazR[i]=v2R[2]+dmz

if (tilta != 0.0 or tiltb != 0.0 or tiltc != 0.0):
    q1R1 = axisangle_to_q(a,np.deg2rad(tilta))
    q1R2 = axisangle_to_q(b,np.deg2rad(tiltb))
    q1R3 = axisangle_to_q(c,np.deg2rad(tiltc))
    [cx1,cy1,cz1] = datamean[0]
    [cx2,cy2,cz2] = datamean[1]
    for j in range(itt):
        i = imap[j] 
        temp = [aax[i]-cx1,aay[i]-cy1,aaz[i]-cz1]
        temp = qv_mult(q1R1,(temp[0],temp[1],temp[2])) 
        temp = qv_mult(q1R2,(temp[0],temp[1],temp[2])) 
        temp = qv_mult(q1R3,(temp[0],temp[1],temp[2])) 
        aax[i]=temp[0]+cx1; aay[i]=temp[1]+cy1; aaz[i]=temp[2]+cz1;        
#
        temp = [aaxR[i]-cx1,aayR[i]-cy1,aazR[i]-cz1]
        temp = qv_mult(q1R1,(temp[0],temp[1],temp[2])) 
        temp = qv_mult(q1R2,(temp[0],temp[1],temp[2])) 
        temp = qv_mult(q1R3,(temp[0],temp[1],temp[2])) 
        aaxR[i]=temp[0]+cx1;aayR[i]=temp[1]+cy1;aazR[i]=temp[2]+cz1;
#
        i = imap[j+itt] 
        temp = [aax[i]-cx2,aay[i]-cy2,aaz[i]-cz2]
        temp = qv_mult(q1R1,(temp[0],temp[1],temp[2])) 
        temp = qv_mult(q1R2,(temp[0],temp[1],temp[2])) 
        temp = qv_mult(q1R3,(temp[0],temp[1],temp[2])) 
        aax[i]=temp[0]+cx2;aay[i]=temp[1]+cy2;aaz[i]=temp[2]+cz2;        
#
        temp = [aaxR[i]-cx2,aayR[i]-cy2,aazR[i]-cz2]
        temp = qv_mult(q1R1,(temp[0],temp[1],temp[2])) 
        temp = qv_mult(q1R2,(temp[0],temp[1],temp[2])) 
        temp = qv_mult(q1R3,(temp[0],temp[1],temp[2])) 
        aaxR[i]=temp[0]+cx2;aayR[i]=temp[1]+cy2;aazR[i]=temp[2]+cz2;
#    print i, xx[i],yy[i],zz[i]
n_list = []
for i in range(l2_num):
    newline = ' '.join(astring[i].split())
    a_list=newline.split(" ")
    m = len(a_list)
    if (m == 1):
        n_list.append([int(a_list[0])])
    elif (m == 2):
        n_list.append([int(a_list[0]),int(a_list[1])])
    elif (m == 3):
        n_list.append([int(a_list[0]),int(a_list[1]),int(a_list[2])])
xnew = []; ynew = []; znew=[]; new_atype = []; new_abt = []; new_astring = []
nstart=0
mynewstring=''
print lenx,leny,lenz
for i in range(0,lenx):
    ai = float(i)
    for j in range(0,leny):
        aj = float(j)
        kstart = 0
        if (rough == 'T'):
            if (i == lenx-2):
                kstart = 1
            elif (i == lenx-1):
                kstart = 2
        for k in range(kstart,lenz):
            ak = float(k)
            nct += l2_num
            flip = 0
            if (i == iflip and j==jflip and k==kflip):
                flip = 1
            elif (iflip == 'all'):
                flip =1
            elif (iflip == 'random'):
                ir = random.randint(1, 2)
                if (ir == 2):
                    flip =  1
            elif (iflip == 'layers'):
                if (i in jflip):
#                    print 'fliping layer ', i
                    flip = 1
            for ii in range(l2_num):
                if (flip == 0):  # a[1]=a[2]=b[2]=0
                    xnew.append(aax[ii] +ai*a[0]+aj*b[0]+ak*c[0])
                    ynew.append(aay[ii] +        aj*b[1]+ak*c[1])
                    znew.append(aaz[ii] +                ak*c[2])
                else:
                    xnew.append(aaxR[ii]+ai*a[0]+aj*b[0]+ak*c[0])
                    ynew.append(aayR[ii]+        aj*b[1]+ak*c[1])
                    znew.append(aazR[ii]+                ak*c[2])
                new_atype.append(atype[ii])
                new_abt.append(abt[ii])
                newstring = ""
                for kk in n_list[ii]:
                    newstring = newstring+"%6s" % (kk+nstart)
#                print ii, newstring
#                raw_input()
                new_astring.append(newstring)
#                if (nsp_flag == 'T' and ii==16 and i == 0): # a dangling H
#                    nspecial += 1
#                    isp = len (xnew)-1
#                    mynewstring="%5s%10.4f%10.4f%10.4f" %  ((nct-l2_num + ii+1),xnew[isp],ynew[isp],znew[isp]) 
#                    nsp.append( mynewstring+ ' 0.1 0.0001'+ '\n' )
#                    nsp.append( (str(nct-l2_num + ii+1)+ '\n') )
            if (nsp_flag == 'T'): # a dangling H
                nspecial = 1
                for ii in range(l2_num):
                    if (abt[ii] == 2):
                        mynewstring= mynewstring +' '+str(ii+1+nct-l2_num) 
            nstart += l2_num
#
# A custom reorientation
#
#il = [0,5];
il = []
for n in il:
    ax=[];ay=[];az=[];imap=[]
    k = l2_num*leny*lenz*n
    for i in range(leny):
        for j in range(lenz):
            for m in range(l2_num):
                ax.append(xnew[k])
                ay.append(ynew[k])
                az.append(znew[k])
                imap.append(k)
                k += 1
#print ay
    ax,ay,az = rebuild_Tx(180.0,0.0,0.,uc,ax,ay,az,l2_num)
    for j in range(len(ax)):
        i=imap[j]
        xnew[i]=ax[j]
        ynew[i]=ay[j]
        znew[i]=az[j]
#print ay
#
#
uc2=[]
uc2.append((ai+1.0)*uc[0]+agap)
uc2.append((aj+1.0)*uc[1]+bgap)
uc2.append((ak+1.0)*uc[2]+cgap)
uc2.append(uc[3])
uc2.append(uc[4])
uc2.append(uc[5])
#
#
#raw_input()
##########################
# Now to save the positions   
##########################
#
write_txyz(txyznewname,nct,new_atype,new_astring,xnew,ynew,znew,new_abt,uc2,header)
#print 'new file:',txyznewname
if (substrate !='none'):
    t2_num,s2_num,xnew,ynew,znew,new_atype,new_astring,new_abt,uc,header = merge_txyz(txyznewname,substrate,ashift,bshift,cshift,special)    
    write_txyz(txyznewname,t2_num+s2_num,new_atype,new_astring,xnew,ynew,znew,new_abt,uc2,header)
#
#
f = open(tnewkey,'w')
f.write('# Force Field Selection \nPARAMETERS        /home/winokur/MolecularTools/ffe/../tinker/params/mm3.prm \n')
f.write('# Crystal Lattice And Periodic Boundary \n\n')
f.write('A-AXIS   '+str(uc2[0])+' \n')
f.write('B-AXIS   '+str(uc2[1])+' \n')
f.write('C-AXIS   '+str(uc2[2])+' \n')
f.write('ALPHA    '+str(uc2[3])+' \n')
f.write('BETA     '+str(uc2[4])+' \n')
f.write('GAMMA    '+str(uc2[5])+' \n')
f.write('\n')
for i in range(nspecial):
#    f.write('RESTRAIN-POSITION '+ str(nsp[i]))
    f.write('#PISYSTEM '+ mynewstring+'\n')
f.close()  # close write file
# fini
