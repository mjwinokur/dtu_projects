# This program recenters the NaT2 monomers above the HOPG and, if desired, creates multilayers
#
import numpy as np
from mw_library import read_txyz
#import math
########################################
#   Main
########################################
#
# This code will only add additional T2 layers to prexisting file of T2+HOPG, or recenter or strip out the HOPG
#
# open a base file of the T2 monomer 
#name='NaT2_2.txyz'  # tinker format with direct structure from cif file
#name='newfile.xyz_'  # tinker format with direct structure from cif file
name='NaT2_min_2x2x1p1_hopg_0bcmisfit.xyz_'
#index=22
index=1
#operation = 'add_layers'
operation = 'strip_hopg'
#operation = 'recenter'
#keyword = 'nohopg'
keyword = 'hopg'  # does the read file have an hopg layer 
readname= name+str(index)
l2_num,aax,aay,aaz,atype,astring,abt,uc,header = read_txyz(readname)
n_add_layers=0
#print n_layers

# Translation matrix for the unit cell
if (uc[0] == 0.0):
    uc[0]=26.0
if (uc[1] == 0.0):
    uc[1] = 5.964
if (uc[2] == 0.0):
    uc[2] = 8.12
uc_ax = 20.552
uc_ay = 0.
uc_az = 0.
uc_bx = 0.
uc_by = uc[1]
uc_bz = 0.
#uc_cx = uc[2]*np.cos(np.pi*uc[4]/180.)
uc_cy = 0.0
#uc_cz = uc[2]*np.sin(np.pi*uc[4]/180.)
uc_cx = 0.0
uc_cy = 0.0
uc_cz = uc[2] #This needs to revised in order to properly deal with monoclinic unit cells using conventional coordinates
# and not use orthorhombic approximates

if (keyword=='hopg'):
# open a 2nd base file of the HOPG substrate 
    with open("p-hopg9x15.xyz", "r") as file:
        text = file.readlines()
        hopg_num = -1
#
#hopgarr = np.empty((86), dtype=object) # Define an empty array
    hopgarr = []
    for line in text:
	line=line.strip();  # get rid of \cr\lf at end of line
        hopg_num += 1
#	print hopg_num, line
	mystring = line
	myhopglist = ' '.join(mystring.split())
	hopgarr.append(myhopglist)
    file.close()  # close read file
else:
    hopg_num=0
n_layers = ((l2_num-hopg_num)/384)-1
#
# Need to recenter NaT2
#
###mylist = hopgarr[hopg_num-1].split(" ")
#hopgx = float(mylist[2])*0.50
#hopgy = float(mylist[3])*0.50
#hopgz = float(mylist[4])*0.50
t2_ct = l2_num - hopg_num
#xcen = ycen = zcen = 0.0
#With the NaT2 monomer 0 to 23 plus 48 to 71 for the same unit
#The from 24 to 47 and the from 72 to 95 form the 2nd unit
#This is true throughout the unit cell  
#There are four pair of chains
mstart = []
mcount = []
nb = 24
nc = 4
for j in range(nc):
    i=j*96
    mstart.append(i)
    mstart.append(i+2*nb)
#        
    mstart.append(i+nb)
    mstart.append(i+3*nb)
    mcount.append(nb)
    mcount.append(nb)
if (n_layers > 0):  # assumes an hopg layer and add those atoms past the hopg
    for k in range(n_layers):
        kct=(k+1)*4*96+hopg_num
        for j in range(nc):
            i = j*96+kct 
	    mstart.append(i)
            mstart.append(i+2*nb)
#        
            mstart.append(i+nb)
            mstart.append(i+3*nb)
	    mcount.append(nb)
	    mcount.append(nb)
#	print 'Layer:',k+1
if (hopg_num>0):
    mstart.append(4*96) # Center the HOPG layer
    mstart.append(4*96+hopg_num/2) # Center the HOPG layer
    mcount.append(hopg_num/2)

k = 0    
ncycle = 2*nc*(n_layers+1)+1
if (hopg_num == 0):
    ncycle=2*nc*(n_layers+1)
for j in range(ncycle):
    xcen = ycen= zcen = 0.0
    i1=mstart[k]
    i3=mstart[k+1]
    nct=mcount[j]
    act=float(nct)
    for i in range(nct):
        i2 = i + i1
        xcen= xcen + aax[i2]
        ycen= ycen + aay[i2]
        zcen= zcen + aaz[i2]
        i4 = i + i3
        xcen= xcen + aax[i4]
        ycen= ycen + aay[i4]
        zcen= zcen + aaz[i4]
    xcen = xcen /act
    ycen = ycen /act
    zcen = zcen /act
    k +=2
#    print uc[0],uc[1],uc[2]
#    print j,act,i1,i3,nct,xcen,ycen,zcen
# So now if the monomer center of is outside the central zone a shift is needed
    xshift= yshift= zshift = 0.0
    if (xcen < 0.0):
        xshift=1.0*uc[0]
    if (xcen > uc[0]):
        xshift= -1.0*uc[0]
    if (ycen < 0.0):
        yshift=1.0*uc[1]
    if (ycen > uc[1]):
        yshift= -1.0*uc[1]
    if (zcen < 0.0):
        zshift=1.0*uc[2]
    if (zcen > uc[2]):
        zshift= -1.0*uc[2]
    if ( abs(xshift) + abs(yshift) + abs(zshift) > 0.0):
	for i in range(nct):
            i2 = i + i1
            aax[i2] = aax[i2] + xshift
            aay[i2] = aay[i2] + yshift
	    aaz[i2] = aaz[i2] + zshift
            i4 = i + i3
            aax[i4] = aax[i4] + xshift
            aay[i4] = aay[i4] + yshift
	    aaz[i4] = aaz[i4] + zshift
#    print  i2,i4,xcen,ycen,zcen,xshift, yshift, zshift

#print l2_num, hopg_num, t2_ct, xcen, ycen, zcen, hopgx, hopgy, hopgz
#
# Now to create multiple NaT2 layers
# 
for i in range(n_add_layers):
    ashift=uc_ax*float(i+1)
    for j in range(t2_ct):
        atype.append(atype[j])
        aax.append(aax[j]+ashift)
	aay.append(aay[j])
	aaz.append(aaz[j])
	newline = ' '.join(astring[j].split())
	mylist = newline.split(" ")
        m = len(mylist)
	temp = 0.0
	newstring = ''
	temp = int(mylist[0])
        newstring = "%6s" % (temp) 

	for n in range(m-1):  #now to put in bonding configuration
            temp = int(mylist[n+1])+l2_num
            tstring = "%6s" % (temp) 
	    newstring=newstring+tstring
#               print(newstring)	
	astring.append(newstring)
    l2_num = l2_num + t2_ct
    print l2_num,t2_ct
    raw_input( )
#
#  Now if need remove hopg
#
if  (operation == 'strip_hopg'):
    l2_num = l2_num - hopg_num
    kct = 96 * 4
    for i in range(n_layers*kct):
	ii = i + kct
	ij = ii + hopg_num
        atype[ii] = atype[ij] 
        aax[ii] = aax[ij] 
	aay[ii] = aay[ij] 
	aaz[ii] = aaz[ij] 
	newline = ' '.join(astring[ij].split())
	mylist = newline.split(" ")
        m = len(mylist)
	temp = 0.0
	newstring = ''
	temp = int(mylist[0])
        newstring = "%6s" % (temp) 

	for n in range(m-1):  #now to put in bonding configuration
            temp = int(mylist[n+1])-hopg_num
            tstring = "%6s" % (temp) 
	    newstring=newstring+tstring
#               print(newstring)	
	astring[ii] =newstring
    uc[0]=133.260000   
    uc[1]=11.928000   
    uc[2]=16.240000

print l2_num,t2_ct
raw_input()
writename= name+str(index+1)
f = open(writename,'w')
write_txyz(writename,l2_num)
#f2 = open(name+'2.xyz','w')

f = open('newfile.key','w')
f.write('# Force Field Selection \nPARAMETERS        /home/winokur/MolecularTools/ffe/../tinker/params/mm3.prm \n')
f.write('# Crystal Lattice And Periodic Boundary \n\n')
if  (operation == 'add_layers'):
    aaxis=str(uc[0]+n_add_layers*(uc_ax)+4.5)
elif  (operation == 'strip_hopg'):
    aaxis=str(uc[0]-9.0)
else:
    aaxis=str(uc[0])
baxis=str(uc[1])
caxis=str(uc[2])
f.write('A-AXIS    '+aaxis+'\n')
f.write('B-AXIS    '+baxis+'\n')    
f.write('C-AXIS    '+caxis+'\n')
f.write('ALPHA             90. \nBETA              90.     \nGAMMA             90. \n')
f.close()  # close write file


