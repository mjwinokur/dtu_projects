#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 13:25:24 2017

@author: winokur
"""
import numpy as np
#import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
from mw_library import iofiles
#from PIL import Image # cv2 is much easier to use
import cv2
import sys
output_file =""
i_file =""
xs=1.0
mask0 = []; mask1 =[]
maska0 = []; maska1 =[]
NaT3_corr_1 = 'false'
NaT3_corr_2 = 'false'
NaT3_corr_3 = 'false'
# Typical use
#python tif_process.py -i /home/winokur/Downloads/nat2_corrected.tiff -o NaT2_minus_bkg.tiff -t "mask 550,672,810,985 610,742,890,1015 integrate 556,701,826 590,735,860 calibrate 0.00195296752 i_file NaT2_ll_data.dat"
#python ../tif_process.py -i sample_1A_101_01364_1408.tiff -o S1A_101_minus_bkg.tiff -t "mask 565,670,820 635,740,890 integrate 575,680,830,465 625,730,880,515 calibrate 0.002038 7 i_file S1A_101_ll_data.dat"
#python ../tif_process.py -i sample_1a_103_00724_0695.tiff -o S1A_103_minus_bkg.tiff -t "mask 565,679,830 645,759,910 integrate 580,694,845,472 654,744,895,522 calibrate 0.002038 13 i_file S1A_101_ll_data.dat"
tiff_file = '/home/winokur/Downloads/nat2_corrected.tiff'
tfile,output_file,text = iofiles(sys.argv[1:])
if (len(tfile) > 0):
    tiff_file=tfile
if (text != ""):
#    print "with special:",text
    alist = text.split()
    j=0
    for i in alist:
        if (alist[j] == 'mask' ):
            mask0 = []; mask1 =[]
            blist = alist[j+1].split(',')
            for k in blist:
                mask0.append(int(k))
            blist = alist[j+2].split(',')
            for k in blist:
                mask1.append(int(k))
        if (alist[j] == 'integrate' ):
            maska0 = []; maska1 =[]
            blist = alist[j+1].split(',')
            for k in blist:
                maska0.append(int(k))
            blist = alist[j+2].split(',')
            for k in blist:
                maska1.append(int(k))
        if (alist[j] == 'i_file' ):
            i_file=alist[j+1]
        if (alist[j] == 'calibrate' ):
            xs=float(alist[j+1])
            ioff=int(alist[j+2])
        if (alist[j] == 'NaT3_corr_1' ):
            NaT3_corr_1 = 'true'
        if (alist[j] == 'NaT3_corr_2' ):
            NaT3_corr_2 = 'true'
        if (alist[j] == 'NaT3_corr_3' ):
            NaT3_corr_3 = 'true'
        j += 1
rc('text', usetex=True)
#
img = cv2.imread(tiff_file,-1) # If <0 this returns a 16 bit image

(nx,ny)=img.shape
x = np.linspace(0, nx, num=nx, endpoint=False)
y = np.linspace(0, ny, num=ny, endpoint=False)
X,Y = np.meshgrid(x,y, copy=False)
# Rather than Z below I need to generate a weighting curve of zeros and ones
#Z = X**2 + Y**2 + np.random.rand(*X.shape)*0.01 # An initial test function
Z = img.astype(float) # convert 16 bit integers to float
W = np.ones((nx,ny))
for i in range(ny):  # Weight regions using only where there is data
    temp = img[:,[i]]
#    temp2 = np.argwhere(temp == 0 )
    temp2 = np.where(temp == 0 ) # find pixels with zero intensity to mask
    W[temp2[0],i]=0.0
W2 = W.astype(float)
spsc = 1.0
# An ad hoc correction for the shadowing of the film
if (NaT3_corr_1 == 'true'):
    print ' Ad hoc NaT3 correction 1'
    spsc = 0.88
    T1 = (230.-60.*Y/800.)/(np.exp((796.-(X-Y*0.082))/4.2)+1.)
    T2 = (1050.+ 70.*Y/650.)/(np.exp((883.-(X-Y*0.086))/9.6)+1.)
    Z -= (T1+T2)*W2
if (NaT3_corr_2 == 'true'):
    print ' Ad hoc NaT3 correction 2'
    spsc = 0.88
    T1 = (72.-15.*Y/800.)/(np.exp((794.-(X-Y*0.082))/4.6)+1.)
    T2 = (95.+ 55.*Y/650.)/(np.exp((883.-(X-Y*0.086))/8.6)+1.)
    Z -= (T1+T2)*W2
if (NaT3_corr_3 == 'true'):
    print ' Ad hoc NaT3 correction 3'
    spsc = 0.88
    T1 = (300.-60.*Y/800.)/(np.exp((794.-(X-Y*0.082))/4.6)+1.)
    T2 = (750.+500.*Y/650.)/(np.exp((883.-(X-Y*0.086))/8.6)+1.)
    Z -= (T1+T2)*W2
#mask0 = [550,672,810,985]
#mask1 = [610,742,890,1015]
for i in range(len(mask0)): # Mask off region to exclude from fitting background
    W[:,mask0[i]:mask1[i]]=0.0
#    print temp2
#    print i
#    plt.plot(x,temp)
#    plt.plot(x,W[:,i])
#    plt.show()
#    raw_input()
#Z = W*Z
X = X.flatten() # convert to 1D arrays
Y = Y.flatten()
Wf = W.flatten()
#
X2 = X**2
Y2 = Y**2
X3 = X**3
Y3 = Y**3
X2Y2 =X2*Y2
XY = X*Y
#A = np.array([X*0+1, X, Y, X2, X2*Y, X2Y2, Y2, X*Y2, XY]).T
A = np.array([X*0+1, X, Y, X2, X2*Y, X2Y2, Y2, X*Y2, XY,X3,Y3,X3*Y3,X*Y3,X3*Y,X2*Y3,X3*Y2]).T
B = Z.flatten()
Aw = A * np.sqrt(Wf[:,np.newaxis]) # add axis for weighting of fit
Bw = B * np.sqrt(Wf)
#Aw = A * Wf[:,np.newaxis]  # doesn't work
#Bw = Wf
#
coeff, r, rank, s = np.linalg.lstsq(Aw, Bw) # fit background
X,Y = np.meshgrid(x,y, copy=False)
X2 = X**2
Y2 = Y**2
X2Y2 =X2*Y2
XY = X*Y
X3 = X**3
Y3 = Y**3
#T2 = (390.-180.*Y/800.)/(np.exp((793.-(X-Y*0.09))/6.)+1.)
#         1. / (exp(-(x-mu)/T) + 1.) 832 step 200 920  
Z2 = (coeff[0]+coeff[1]*X+coeff[2]*Y+coeff[3]*X2+coeff[4]*X2*Y+coeff[5]*X2Y2+coeff[6]*Y2+
     coeff[7]*X*Y2+coeff[8]*XY)
Z2 += (coeff[9]*X3+coeff[10]*Y3+coeff[11]*X3*Y3+coeff[12]*X*Y3+coeff[13]*X3*Y+coeff[14]*X2*Y3+coeff[15]*X3*Y2)
#Z2 += (coeff[16]*T2)
Z3 = (Z-spsc*Z2*W2)  # some simple renormalizations


Z3[Z3 <0]=0.1  #  To limit display artifacts
show_mask = 'true'
if (show_mask == 'true'):
#    plt.subplot(121),plt.imshow(np.sqrt(Z), cmap = 'terrain')  # data
#    plt.title('Original Image ($\sqrt I$) '), plt.xticks([]), plt.yticks([])
#    plt.subplot(121),plt.imshow((Z2*W2), cmap = 'terrain')
#    plt.title('Calc. background profile '), plt.xticks([]), plt.yticks([])
#    plt.subplot(121),plt.imshow((T2), cmap = 'terrain')
    plt.subplot(121),plt.imshow((Z*W), cmap = 'terrain')
    plt.title('Original Image and Mask ') #, plt.xticks([]), plt.yticks([])
    plt.subplot(122),plt.imshow(W, cmap = 'terrain')
    
#    plt.subplot(122),plt.imshow(np.sqrt(Z3+1.)-1., cmap = 'terrain')
    #plt.subplot(122),plt.imshow(np.sqrt(Z3)-4., cmap = 'terrain')
    plt.title('Masked regions') #, plt.xticks([]), plt.yticks([])
    plt.show()

#plt.subplot(121),plt.imshow(np.sqrt(Z), cmap = 'terrain')  # data
plt.subplot(121),plt.imshow(Z, cmap = 'terrain')  # data
plt.title('Original Image'), plt.xticks([]), plt.yticks([])
#plt.subplot(122),plt.imshow(W, cmap = 'terrain')
#plt.subplot(122),plt.imshow((Z2*W2), cmap = 'terrain')

plt.subplot(122),plt.imshow(np.sqrt(Z3+1.)-1., cmap = 'terrain')
#plt.subplot(122),plt.imshow(np.sqrt(Z3)-4., cmap = 'terrain')
plt.title('Modified Spectrum ($\sqrt I$)') , plt.xticks([]), plt.yticks([])
plt.show()

fft='false'
if (fft == 'true'):
    f = np.fft.fft2(Z3)
    fshift = np.fft.fftshift(f)
    magnitude_spectrum = 20*np.log(np.abs(fshift))
    #plt.subplot(121),plt.imshow(img, cmap = 'gray')
    plt.subplot(121),plt.imshow(Z3, cmap = 'terrain')
    plt.title('Input Image'), plt.xticks([]), plt.yticks([])
    plt.subplot(122),plt.imshow(magnitude_spectrum, cmap = 'terrain')
    plt.title('Magnitude Spectrum'), plt.xticks([]), plt.yticks([])
    plt.show()
#

# Need to make this more IO oriented
#maska0 = [556,701,826]
#maska1 = [590,735,860]
maska2 = [10/35.,10./35.,10./35.,10./35.]
if (len(maska0) >0):
#    ioff = 7
    xoff = len(x)-ioff
#xp = np.linspace(0, xoff, num=xoff, endpoint=False)
    xp = np.zeros(xoff)
#    xs=1.744/893.
#    xs=0.00195296752
    datas = []
    for i in range(xoff):
        xp[i]=float(i)*xs
        newstring = "%7.4f" % xp[i]
        datas.append(newstring)
    for j in range(len(maska0)):
        scan =[]
        for i in range (len(Z3[j])):
            scan.append(Z3[i,maska0[j]:maska1[j]].sum())
        scan0=np.asanyarray(scan[ioff:])/float(maska1[j]-maska0[j]+1)
        plt.plot(xp,scan0*maska2[j])
#        plt.plot(xp,np.sqrt(scan0)*maska2[j])
        for i in range(xoff):
            newstring = "%12.4f" % scan0[i]
            datas[i] += ' '+newstring
    plt.title('Bragg rod integrated intensities') #, plt.xticks([]), plt.yticks([])
    if (i_file != ""):
        print "Writing Bragg rod scans: ",i_file
        f = open(i_file,'w')
        for i in range(xoff):
            f.write(datas[i]+'\n')
        f.close()
#    plt.xticks([0])
#    plt.yticks([])
    plt.xlim([0.,max(xp)])
    plt.xlabel('Wavevector')
    plt.ylabel('Intensity')
#    plt.axes()
    plt.show()

#
if (output_file != ''):
    Z2 = np.uint16(Z3) # convert to unsigned 16 bit format
#Z = Z3.astype(int)
    print "Writing tiff output file to:", output_file
    cv2.imwrite(output_file,Z2)
    Z2 = np.uint16(Z3) # convert to unsigned 16 bit format
#Z = Z3.astype(int)
    Z2 = np.sqrt(Z3+1.)-1. # convert to unsigned 16 bit format
    Z2 = np.uint16(Z2) # convert to unsigned 16 bit format
    output_file = 'sqrt_'+output_file
    print "Writing tiff output file to:", 'sqrt_'+output_file
    cv2.imwrite(output_file,Z2)
