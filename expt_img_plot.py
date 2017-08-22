#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 14:26:16 2017

@author: winokur
"""
from __future__ import print_function
import numpy as np
#import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
from mw_library import iofiles
#from PIL import Image # cv2 is much easier to use
import cv2
import sys
# from scipy.interpolate import griddata
#import sys
# because of the small number of lattice repeats a Lorentzian is used to 
# Now to plot the experimental tiff files
#
image_dir = '/home/winokur/dtu_projects/NaT3_sample/'
tiff_file = 'sqrt_NaT3_02272_2744_minus_bkg.tiff'
image_source = image_dir+tiff_file
img = cv2.imread(image_source,-1) # If <0 this returns a 16 bit image
my_color='terrain'
#interp='bilinear'
interp='spline36'
xmin = 0
xmax = 1024
ymin = 0
ymax = 1024

(nx,ny)=img.shape
x = np.linspace(0, nx, num=nx, endpoint=False)
y = np.linspace(0, ny, num=ny, endpoint=False)
X,Y = np.meshgrid(x,y, copy=False)
Z = img.astype(float) # convert 16 bit integers to float

im = plt.imshow(Z, interpolation=interp, origin='lower', cmap = my_color, extent=(xmin,xmax, ymin,ymax))

raw_input()
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

levels = np.arange(300.,1700., 400.)
custom_legend='''NaT3 ($\sqrt I$) '''
#custom_legend=''
#                 colors=('r', 'green', 'blue', (1, 1, 0), '#afeeee', '0.5')
#CS = plt.contour(ZJ, levels,colors=('#afcccc', (1, 1, 0), '#afeeee', '0.5'),
#                 origin='lower',
#                 linewidths=0.5,
#                 extent=(xmin,xmax,ymin,ymax))
# Thicken the zero contour (at 100.)
#zc = CS.collections[0]  # Choose the particular contour for special sizing

#plt.setp(zc, linewidth=0.5)

#plt.clabel(CS, levels[1::3],  # label every third level
#           inline=1,fmt='%1.0f',fontsize=5)
# make a colorbar for the contour lines
#CB = plt.colorbar(CS, shrink=1.0, extend='both')
if (custom_legend != ''):
    plt.title(custom_legend)
else:
    plt.title(r'NaT2 or NaT3 test calculation ($\sqrt I$~)')
#plt.hot()  # Now change the colormap for the contour lines and colorbar
plt.xlabel(r'$q_{\rm xy}$ ($\mbox{\AA}^{-1}$) ',fontsize=14)
plt.ylabel(r'$q_{\rm z}$ ($\mbox{\AA}^{-1}$) ',fontsize=14)
plt.flag()
# We can still add a colorbar for the image, too.
CBI = plt.colorbar(im, orientation='vertical', shrink=0.7)
# This makes the original colorbar look a bit out of place,
# so let's improve its position.
l, b, w, h = plt.gca().get_position().bounds
#ll, bb, ww, hh = CB.ax.get_position().bounds
#CB.ax.set_position([ll, b + 0.2*h, 1.2*ww, h*0.3])
#for i in range(len(label_ind)):
#    plt.text(label_x[i],label_y[i], label_ind[i])
#texts = []
#for x, y, text in zip(label_x, label_y, label_ind):
#    texts.append(plt.text(x, y, text, color='y',fontsize=11))
#adjust_text(texts, precision=0,force_text=0.01, size=7,autoalign='y',only_move={'points':'x', 'text':'x'},color='yellow')
