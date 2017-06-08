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
#from PIL import Image # cv2 is much easier to use
import cv2
rc('text', usetex=True)
#
int_type='log'
int_type='sqrt'
colortype='gray'
#colortype='terrain'
itype = 'cat'
itype = 'circle'
if (itype == 'cat'):
    cat_file = 'cat.png'
    img = cv2.imread(cat_file,0) # If <0 this returns a 16 bit image
# Make an image
(nx,ny)=(512,512)
x = np.linspace(0, nx, num=nx, endpoint=False)
y = np.linspace(0, ny, num=ny, endpoint=False)
X,Y = np.meshgrid(x,y, copy=False)
# Rather than Z below I need to generate a weighting curve of zeros and ones
#Z = X**2 + Y**2 + np.random.rand(*X.shape)*0.01 # An initial test function
ii=2
nm=2
(ixr,iyr)=(nm*4*64,nm*4*64)
xr=nx/ixr
yr=ny/iyr
Za= np.zeros([ixr,iyr])+4.
#Z = np.zeros((nx,ny))+0.0
if (itype == 'circle'):
    l1=40
    l2=60
    radius2=l2*l2
    radius1=l1*l1
    for i in range(ixr/2-l2-1,ixr/2+l2+1):
        for j in range(iyr/2-l2-1,iyr/2+l2+1):
            rcal=(i-ixr/2)**2.+(j-iyr/2)**2.
            if (rcal < radius1):
                Za[i,j]=4.0
            elif (rcal < radius2):
                Za[i,j]=8.
elif (itype == 'cat'):
    ixoff=nx/2-2*nx/7/nm
    iyoff=ny/2-2*ny/5/nm
    (ix,iy) = (img.shape)
    catmap = np.ones([ix,iy])
    for i in range(ix):
        temp = img[:,[i]]
#    temp2 = np.argwhere(temp == 0 )
        temp2 = np.where(temp == 255 ) # find pixels with zero intensity to mask
        catmap[temp2[0],i]=0.0
        temp2 = np.where(temp != 255 ) # find pixels with zero intensity to mask        
        Za[temp2[0]+iyoff,i+ixoff]=1.
    
    
ix = np.linspace(0, nx, num=8, endpoint=False)
iy = np.linspace(0, ny, num=8, endpoint=False)
Zp = np.tile(Za,xr)
Zpp = Zp
for i in range(1,yr):
    Zpp = np.vstack((Zpp,Zp))
#Z3 = (Z-Z2*W2)  # +16.) some simple renormalizations
#Z3[Z3 <0]=0.1
###plt.subplot(121),plt.imshow(Zpp, cmap = 'terrain')  # data
#plt.title('Original Image'), plt.xticks([]), plt.yticks([])
#plt.title('Modified Spectrum ($\sqrt I$)'), plt.xticks([]), plt.yticks([])
###plt.show()

f = np.fft.fft2(Zpp)
#f = np.fft.rfft2(Zpp)
fshift = np.fft.fftshift(f)
fphase = np.angle(f)

if (int_type == 'log'):
    magnitude_spectrum = 10.*np.log(np.abs(fshift))
    title = 'FFT Intensity (log)'
elif (int_type == 'sqrt'):
    magnitude_spectrum = 10.*np.sqrt(np.abs(fshift))
    title = 'FFT Intensity (sqrt)'

#plt.subplot(121),plt.imshow(img, cmap = 'gray')
plt.subplot(121),plt.imshow(Zpp, cmap = colortype)
plt.title('Input Image'), plt.xticks([]), plt.yticks([])
if (itype == 'circle'):
    for i in range(nx/2-4,nx/2+4):
        for j in range(ny/2-4,ny/2+4):
            magnitude_spectrum[i,j]=0.

plt.subplot(122),plt.imshow(magnitude_spectrum, cmap = colortype, interpolation='gaussian')
plt.title('FFT Intensity (log)'), plt.xticks([]), plt.yticks([])
#plt.subplot(122),plt.imshow(magnitude_spectrum, cmap = 'terrain', interpolation='bicubic')
#plt.title('Magnitude Spectrum'), plt.xticks([]), plt.yticks([])
plt.show()
plt.subplot(121),plt.imshow(Zpp, cmap = colortype)
plt.title('Input Image'), plt.xticks([]), plt.yticks([])
plt.subplot(122),plt.imshow(fphase*10, cmap = colortype, interpolation='bicubic')
plt.title('FFT phase'), plt.xticks([]), plt.yticks([])
#plt.subplot(122),plt.imshow(magnitude_spectrum, cmap = 'terrain', interpolation='bicubic')
#plt.title('Magnitude Spectrum'), plt.xticks([]), plt.yticks([])
plt.show()
Z3 = np.uint16(Zpp) # convert to unsigned 16 bit format
#Z = Z3.astype(int)
#cv2.imwrite('image.tiff',Z3)
