#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 10:35:01 2017

@author: winokur
"""

import matplotlib
import numpy as np
import math
#from scipy.special import i0
#import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
from matplotlib import rc
from mw_library import read_txyz_info_uc
from mw_library import abc_to_rlp
from mw_library import read_hklI_file
from adjustText import adjust_text
from mw_library import iofiles
from mw_library import read_data_file  # Experimental data
#from scipy.interpolate import griddata
import sys
tinkfile,hklfile,text = iofiles(sys.argv[1:])
store = "false"
# Save to an output list and check for a custom input list
sigmax0=0.006; sigmay0=0.009; sigma_arc=0.015 # in radians
alpha=0.; beta=0.; gamma = 0.0
lldata = 'none'
hc = []; kc =[]; lc=[]
if (text != ""):
    print "with special:",text
    alist = text.split()
    print alist
    j=0
    for i in alist:
#        print j,i
        if (i == 'hc'):
            blist = alist[j+1].split(',')
            for k in blist:
                hc.append(int(k))
        if (i == 'kc'):
            blist = alist[j+1].split(',')
            for k in blist:
                kc.append(int(k))
        if (i == 'lc'):
            blist = alist[j+1].split(',')
            for k in blist:
                lc.append(int(k))
        if (i == 'store'):
            store = 'true'
            sfile = alist[j+1]
        elif (i == 'sigma_arc'):
            sigma_arc = float(alist[j+1])
        elif (i == 'sigmax'):
            sigmax0 = float(alist[j+1])
        elif (i == 'sigmay'):
            sigmay0 = float(alist[j+1])
        elif (alist[j] == 'alpha'): # An arbitrary space
            alpha = float(alist[j+1])
            print 'alpha: ',alpha
        elif (alist[j] == 'beta'): # An arbitrary space
            beta = float(alist[j+1])
            print 'beta: ',beta
        elif (alist[j] == 'gamma'): # An arbitrary space
            gamma = float(alist[j+1])
            print 'gamma: ',gamma
        elif (alist[j] == 'data'): # An arbitrary space
            lldata = alist[j+1]
            print 'Expt. data file: ',lldata
        j += 1

    if (len(hc) != len(kc) or len(hc) != len(lc) or len(kc) != len(lc)):
        print hc, lc, kc
        print 'hc kc lc must contain an identical number of elements'
        raw_input()
rc('text', usetex=True)
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
# Load in sf calculation set
#
if (tinkfile == '' or hklfile == ''):
    if (tinkfile == '' and hklfile == ''):
# Decode a tinker file
##tinkfile = '/home/winokur/Documents/gixrd/NaT2mw.txyz'
##hklfile = '/home/winokur/Documents/gixrd/NaT2_mw.hkl'
        tinkfile = '/home/winokur/dtu_projects/ztest.100'  # At 300 K
        hklfile = '/home/winokur/dtu_projects/ztest.hkl'
        tinkfile = '/home/winokur/dtu_projects/ztest.300'  # At 400 K
        hklfile = '/home/winokur/dtu_projects/ztest.hkl_300'
        tinkfile = '/home/winokur/dtu_projects/ztest.500'  # At 500 K
        hklfile = '/home/winokur/dtu_projects/ztest.hkl_500'
        tinkfile = '/home/winokur/dtu_projects/ztest.700'  # At 550 K
        hklfile = '/home/winokur/dtu_projects/ztest.hkl_700'
        tinkfile = '/home/winokur/dtu_projects/ztest.900'  # At 600 K
        hklfile = '/home/winokur/dtu_projects/ztest.hkl_900'
        tinkfile = '/home/winokur/dtu_projects/ztest.1100'  # At 700 K
        hklfile = '/home/winokur/dtu_projects/ztest.hkl_1100'
        tinkfile = '/home/winokur/dtu_projects/ztest.1300'  # At 800 K
        hklfile = '/home/winokur/dtu_projects/ztest.hkl_1300'
        tinkfile = '/home/winokur/dtu_projects/ztest.1500'  # At 900 K
        hklfile = '/home/winokur/dtu_projects/ztest.hkl_1500'
        tinkfile = '/home/winokur/dtu_projects/z2test.200'  # At 900 K
        hklfile = '/home/winokur/dtu_projects/z2test.hkl_200'
        tinkfile = '/home/winokur/dtu_projects/z2test.xyz_3'  # At 0 K
        hklfile = '/home/winokur/dtu_projects/z2test.hkl_xyz_3'
        tinkfile = '/home/winokur/dtu_projects/z2test.xyz_5'  # At 0 K new unit cell
        hklfile = '/home/winokur/dtu_projects/z2test.hkl_xyz_5'
#tinkfile = '/home/winokur/dtu_projects/ztest.xyz' # Minimumized at 0 K
#hklfile = '/home/winokur/dtu_projects/ztest0.hkl'
#tinkfile = '/home/winokur/dtu_projects/wtest.xyz_14' # Minimumized at 0 K
#hklfile = '/home/winokur/dtu_projects/wtest14.hkl'        
    else:
        print('Both the input and the pseudo output files are needed!')
        raw_input()
#
indfile='index_out.txt'
indfile_input=''
#
uc = read_txyz_info_uc(tinkfile)
if (len(uc)== 0):
    raw_input('Stop, this tinker file requires unit cell information')
if (alpha != 0.):
    uc[3]=alpha
if (beta != 0.):
    uc[4]=beta
if (gamma != 0.):
    uc[5]=gamma
# From a,b,c, alpha, beta,gamma get the x,y,z components
astar,bstar,cstar = abc_to_rlp(uc)
astarhat=astar/np.sqrt(np.dot(astar,astar))
#bstarhat=bstar/np.sqrt(np.dot(bstar,bstar))
#cstarhat=cstar/np.sqrt(np.dot(cstar,cstar))
# deduce k_para and k_perp 
hkl,I = read_hklI_file(hklfile)
l_num = len(hkl)
Gper = []; Gpar = []; I_k = []
#xdspmax=2.; xdspmin=0.0 # display limits
#ydspmax=2.; yspdmin=0.0 #
xmax=2.; xmin=0.0 # limits over which the calculation must be done
ymax=2.; ymin=0.0 #
#
label_x = []; label_y = []; label_ind = [];ind_tex = []
j = 30
ind = np.argsort(I)[-j:]
label_xoff = (xmax-xmin)/32.
label_yoff = (ymax-ymin)/100.
for i in ind:
    [hh,kk,ll]=hkl[i]
    Ghkl = hh*astar+kk*bstar+ll*cstar
# a* must be perpendicular to b and c or parallel to b X c
    Gperp = abs(np.dot(Ghkl,astarhat))
    Gmag2 = np.dot(Ghkl,Ghkl)
    Gpara= np.sqrt(abs(Gmag2-Gperp*Gperp))  # can be just less than zero
    ih,ih2 = divmod(hh,3);ik,ik2 = divmod(kk,3); il,il2 = divmod(ll,3);
    ih = int(ih);ik = int(ik); il = int(il);
#    ih = int(hh);ik = int(kk); il = int(ll);
    if (Gperp < ymax and Gpara < xmax):
#        print '%4s%3s%3s%3s%1s%8.2f%1s%6.3f%6.3f' % (i+2,ih,ik,il,' ',I[i],' ',Gpara,Gperp)
        # x and y are transposed on plotting
##        Gpara += label_xoff
##       Gperp -= label_yoff
#            print 'here'
        label_x.append(Gpara)
        label_y.append(Gperp)
        ith = ''
        if (ih2 != 0.):
            ith = '.'+str(int(10.*ih2/3.)) 
        if (ih < 0):
            hlabel='\overline{'+str(-ih)+ith+'}'
        else:
            hlabel= str(ih)+ith
        itk = ''
        if (ik2 != 0.):
            itk = '.'+str(int(10.*ik2/3.)) 
        if (ik < 0):
            klabel='\overline{'+str(-ik)+itk+'}'
        else:
            klabel= str(ik)+itk
        itl = ''
        if (il2 != 0.):
            itl = '.'+str(int(10.*il2/3.)) 
        if (il < 0):
            llabel='\overline{'+str(-il)+itl+'}'
        else:
            llabel= str(il)+itl
        label_ind.append('$ ('+hlabel+'\,'+klabel+'\,'+llabel+') $')
        if (ith == '' and itk == '' and itl == ''):
            ind_tex.append( '%3s%3s%3s%1s%8.2f%1s%6.3f%6.3f' % (ih,ik,il,' ',I[i],' ',Gpara,Gperp))
        else:
            temp ='%3s%3s%3s%1s%8.2f%1s%6.3f%6.3f' % (ih,ik,il,' ',I[i],' ',Gpara,Gperp)
            temp = temp+' *** '+str(ih)+ith+' '+str(ik)+itk+' '+str(il)+itl
            ind_tex.append(temp)
f = open(indfile,'w')
for i in ind_tex:
    f.write(i+' \n')
f.close()  # close write file

Narc_pts = int(6.*max([sigmax0,sigmay0,sigma_arc])/0.014)
#Narc_pts = 6 # Total number is 2*Narc_pts-1
Imaxdiv500 = max(I)/500.
sigmax=[];sigmay=[];
if (Narc_pts > 1):
    Iarc_ang = np.linspace(-3.*sigma_arc,3.*sigma_arc,num=(Narc_pts*2-1))
    Iarc_int = 1./(sigma_arc*np.sqrt(2*np.pi))*np.exp(-Iarc_ang**2./(2.*sigma_arc*sigma_arc) )
#    Iarc_int = np.exp(sigma_arc*np.cos(Iarc_ang))/(2*np.pi*i0(sigma_arc))/area
    Iarc_int = Iarc_int/np.sum(Iarc_int)  # Coarse discrete distributions break the normalizations
###################  NOTE:  Need to check to see if the bivariate normal breaks as well. ######################
    Narc_pts = Narc_pts*2-1 # Total number is 2*Narc_pts-1
#    plt.plot(Iarc_ang,Iarc_int, linewidth=2, color='r')
#    plt.show()
for i in range(l_num):
    [hh,kk,ll]=hkl[i]
    Ghkl = hh*astar+kk*bstar+ll*cstar
    Gperp = abs(np.dot(Ghkl,astarhat))
    Gmag2 = np.dot(Ghkl,Ghkl)
    Gpara= np.sqrt(abs(Gmag2-Gperp*Gperp))  # can be just less than zero       
    scalex = 1.0
    scaley = 1.0
    if (Gpara <0.05):  # narrower widths for reflections along the Gperp direction, should be a more general function
        scaley =0.4
    Gparmax = xmax+0.9*4.*sigmax0*scalex
    Gpermax = ymax+0.9*4.*sigmay0*scaley
#    print i, hkl[i],Gpara,Gperp
#    if (20*int(i/20) == i):
#     raw_input()
    if (Narc_pts == 1 or I[i] < Imaxdiv500):
        if (Gpara < Gparmax and Gperp < Gpermax):
            sigmax.append(sigmax0*scalex)
            sigmay.append(sigmay0*scaley)
            Gpar.append(Gpara)        
            Gper.append(Gperp)
            I_k.append(I[i])
#          
    else:
        Gradial = math.sqrt(Gpara*Gpara+Gperp*Gperp)
        Gangle= math.atan2(Gperp,Gpara)
        for j in range(Narc_pts):        
            temp_ang = Gangle+Iarc_ang[j]
            Gpara=Gradial*np.cos(temp_ang)
            Gperp=Gradial*np.sin(temp_ang)
            if (Gpara < Gparmax and Gperp < Gpermax and Gpara >= -0.000001 and Gperp >= -0.000001):
                sigmax.append(sigmax0*scalex)
                sigmay.append(sigmay0*scaley)
                Gpar.append(Gpara)        
                Gper.append(Gperp)
                I_k.append(I[i]*Iarc_int[j])
# d = 2.*pi/np.sqrt(Gmag2)
# interpolate on k-spacing
# for a particular parallel and perp resolution one needs to 
# convolute 2D Gaussian
# plot
delta = 0.0125; 
xmax2 = xmax-delta; ymax2 = ymax-delta;
#x = np.arange(-3.0, 3.0, delta); y = np.arange(-2.0, 2.0, delta)
# linspace is said to give more consistent results  
Nx=int((xmax-xmin)/delta)  # actually x-y
Ny=int((ymax-ymin)/delta)  # actually z
Z = np.zeros(shape=(Nx,Ny))
xtot = np.linspace(xmin,xmin,Nx+1)
ytot = np.linspace(ymin,ymin,Ny+1)
fluct = 4.0  # how far to span the bivariate_normal function
#Zsub = np.array((Nsx, Nsy))
#
# As it now stands the k_perp (y) border is not properly dealt with because
# any scattering below k_perp = 0.0 is not reflect in the maps.  In principle
# we need to offset the arrays and then map the result onto a display frame
# mirroring the data below the zero.
#
for k in range (len(Gpar)):
    xval=Gpar[k]
    yval=Gper[k]
    tsigmax=2.0*fluct*sigmax[k]
    tsigmay=2.0*fluct*sigmay[k]
    Nsx=int(tsigmax/delta)
    Nsx2=(Nsx+1)/2
    Nsy=int(tsigmay/delta)
    Nsy2=(Nsy+1)/2
    xs_min = xval-fluct*sigmax[k]
    xs_max = xval+fluct*sigmax[k]
    ys_min = yval-fluct*sigmay[k]
    ys_max = yval+fluct*sigmay[k]
    if ((xs_min <  xmin) or (xs_max >  xmax) or (ys_min <  ymin) or (ys_max >  ymax)):
        Nxoff=Nsx
        Nyoff=Nsy
        if (xs_min <  xmin):
            Nxoff=Nxoff-int((xmin-xs_min)/delta)
            xs = np.linspace(xmin,xmin+delta*Nxoff,num=Nxoff, endpoint=False)
            Nsubx=0
        else:
            Nsubx=int((xval-xmin)/delta)-Nsx2
            xb = Nsubx*delta
            xs = np.linspace(xb,xb+tsigmax,num=Nsx, endpoint=False)            
        if (ys_min <  ymin):
            Nyoff=Nyoff-int((ymin-ys_min)/delta)
            ys = np.linspace(ymin,ymin+delta*Nyoff,num=Nyoff, endpoint=False)
            Nsuby=0
        else:
            Nsuby=int((yval-ymin)/delta)-Nsy2
            yb = Nsuby*delta
            ys = np.linspace(yb,yb+tsigmay,num=Nsy, endpoint=False)
        if (xs_max >  xmax2):
            Nsubx=int((xval-xmin)/delta)-Nsx2
            xb = Nsubx*delta
            Nxoff=Nx-Nsubx
            xs = np.linspace(xb,xmax2,num=Nxoff, endpoint=False)
#            Nxoff=Nsx-int((xs_max-xmax2)/delta)
#            xs = np.linspace(xb,xb+delta*Nxoff,num=Nxoff, endpoint=False)
        elif (xs_min >= xmin):
            Nsubx=int((xval-xmin)/delta)-Nsx2
            xb = Nsubx*delta
            xs = np.linspace(xb,xb+tsigmax,num=Nsx, endpoint=False)
        if (ys_max >  ymax):
            Nsuby=int((yval-ymin)/delta)-Nsy2
            yb = Nsuby*delta
            Nyoff=Ny-Nsuby
            ys = np.linspace(yb,ymax2,num=Nyoff, endpoint=False)
#            Nyoff=Nsy-int((ys_max-ymax2)/delta)
#            ys = np.linspace(yb,yb+delta*Nyoff,num=Nyoff, endpoint=False)
        elif (ys_min >= ymin):
            Nsuby=int((yval-ymin)/delta)-Nsy2
            yb = Nsuby*delta
            ys = np.linspace(yb,yb+tsigmay,num=Nsy, endpoint=False)
        XS, YS = np.meshgrid(xs, ys)
        Zsub = mlab.bivariate_normal(XS, YS, sigmax[k], sigmay[k], xval, yval)
        Nxmax=Nxoff
        if ((Nxmax+Nsubx) == Nx):  # A backhanded fix to an edge issue
            Nxmax -= 1
        Nymax=Nyoff
        if ((Nymax+Nsuby) == Ny):
            Nymax -= 1
        for i in range(Nxmax):
            Nxs=i+Nsubx
            for j in range(Nymax):
                Nys=j+Nsuby
                Z[Nxs,Nys] = Z[Nxs,Nys]+Zsub[j,i]*I_k[k]
# This bit can probably be integrated with the above bit but, for the moment it stays
    else:
        Nsubx=int((xval-xmin)/delta)-Nsx2
        xb = Nsubx*delta
        Nsuby=int((yval-ymin)/delta)-Nsy2
        yb = Nsuby*delta
        xs = np.linspace(xb,xb+tsigmax,num=Nsx, endpoint=False)
        ys = np.linspace(yb,yb+tsigmay,num=Nsy, endpoint=False)
        XS, YS = np.meshgrid(xs, ys)
        Zsub = mlab.bivariate_normal(XS, YS, sigmax[k], sigmay[k], xval, yval)
        for i in range(Nsx):
            Nxs=i+Nsubx
            for j in range(Nsy):
                Nys=j+Nsuby
#                print i,j,Nxs,Nys,Nx,Ny
#                raw_input()
                Z[Nxs,Nys] = Z[Nxs,Nys]+Zsub[j,i]*I_k[k]
    if (500*int(k/500) == k): 
        print 'At reflection',k,'/',len(Gpar)
#Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
# define grid.
xi = np.linspace(xmin,xmax,Nx,endpoint=False )
yi = np.linspace(ymin,ymax,Ny,endpoint=False)
XI, YI = np.meshgrid(xi, yi)
ZI=np.sqrt(Z)
ZI=np.transpose(ZI)
#    zi = griddata((XS, YS), Z, (xi[None,:], yi[:,None]), method='cubic')
# Create a simple contour plot with labels using default colors.  The
# inline argument to clabel will control whether the labels are draw
# over the line segments of the contour, removing the lines beneath
# the label
"""
plt.figure()
CS = plt.contour(XI, YI, ZI)
#CS = plt.contour(XS, YS, Zsub)
plt.clabel(CS, inline=1, fontsize=6)
plt.title('Simplest default with labels')
"""
# Or you can use a colormap to specify the colors; the default
# colormap will be used for the contour lines
plt.figure(1,figsize=(8,8)) # default seem to be centimeters
# [None, 'none', 'nearest', 'bilinear', 'bicubic', 'spline16',
#           'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric',
#           'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']
#Possible values are: Accent, Accent_r, Blues, Blues_r, BrBG, BrBG_r, BuGn, BuGn_r, BuPu, BuPu_r, CMRmap, CMRmap_r, 
#Dark2, Dark2_r, GnBu, GnBu_r, Greens, Greens_r, Greys, Greys_r, OrRd, OrRd_r, Oranges, Oranges_r, PRGn, PRGn_r, 
#Paired, Paired_r, Pastel1, Pastel1_r, Pastel2, Pastel2_r, PiYG, PiYG_r, PuBu, PuBuGn, PuBuGn_r, PuBu_r, PuOr, PuOr_r, 
#PuRd, PuRd_r, Purples, Purples_r, RdBu, RdBu_r, RdGy, RdGy_r, RdPu, RdPu_r, RdYlBu, RdYlBu_r, RdYlGn, RdYlGn_r, Reds, 
#Reds_r, Set1, Set1_r, Set2, Set2_r, Set3, Set3_r, Spectral, Spectral_r, Vega10, Vega10_r, Vega20, Vega20_r, Vega20b, 
#Vega20b_r, Vega20c, Vega20c_r, Wistia, Wistia_r, YlGn, YlGnBu, YlGnBu_r, YlGn_r, YlOrBr, YlOrBr_r, YlOrRd, YlOrRd_r, 
#afmhot, afmhot_r, autumn, autumn_r, binary, binary_r, bone, bone_r, brg, brg_r, bwr, bwr_r, cool, cool_r, coolwarm, 
#coolwarm_r, copper, copper_r, cubehelix, cubehelix_r, flag, flag_r, gist_earth, gist_earth_r, gist_gray, gist_gray_r, gist_heat, 
#gist_heat_r, gist_ncar, gist_ncar_r, gist_rainbow, gist_rainbow_r, gist_stern, gist_stern_r, gist_yarg, gist_yarg_r, 
#gnuplot, gnuplot2, gnuplot2_r, gnuplot_r, gray, gray_r, hot, hot_r, hsv, hsv_r, inferno, inferno_r, jet, jet_r, magma, magma_r, 
#nipy_spectral, nipy_spectral_r, ocean, ocean_r, pink, pink_r, plasma, plasma_r, prism, prism_r, rainbow, rainbow_r, seismic, 
#seismic_r, spectral, spectral_r, spring, spring_r, summer, summer_r, terrain, terrain_r, viridis, viridis_r, winter, winter_r
#my_color='Greys_r'
my_color='terrain'
#interp='bilinear'
interp='spline36'
#ax = [plt.subplot(gs[0]),plt.subplot(gs[1]),plt.subplot(gs[2])]
h, w = ZI.shape
#print h,w 125,125
#gs = gridspec.GridSpec(2, 2,width_ratios=[w,w*.2], height_ratios=[h,h*.2])
#im = [plt.subplot(gs[0]),plt.subplot(gs[1])]
#ax[0].imshow(data, cmap='gray', extent = bounds, origin='lower')
#ax[1].plot(X[h/2,:],data[h/2,:],'.',X[h/2,:],data[h/2,:])

im = plt.imshow(ZI, interpolation=interp, origin='lower',
                cmap=my_color, extent=(xmin,xmax, ymin,ymax))
levels = np.arange(100.,1700., 400.)
#                 colors=('r', 'green', 'blue', (1, 1, 0), '#afeeee', '0.5')
CS = plt.contour(ZI, levels,colors=('#afcccc', (1, 1, 0), '#afeeee', '0.5'),
                 origin='lower',
                 linewidths=0.5,
                 extent=(xmin,xmax,ymin,ymax))
# Thicken the zero contour (at 100.)
zc = CS.collections[0]  # Choose the particular contour for special sizing
plt.setp(zc, linewidth=0.5)

plt.clabel(CS, levels[1::3],  # label every third level
           inline=1,
           fmt='%1.0f',
           fontsize=5)
# make a colorbar for the contour lines
CB = plt.colorbar(CS, shrink=1.0, extend='both')
plt.title(r'NaT2 or NaT3 test calculation ($\sqrt I$~)')
#plt.hot()  # Now change the colormap for the contour lines and colorbar
plt.xlabel(r'k$_{xy}$ ($\mbox{\AA}^{-1}$) ')
plt.ylabel(r'k$_z$ ($\mbox{\AA}^{-1}$) ')
plt.flag()
# We can still add a colorbar for the image, too.
CBI = plt.colorbar(im, orientation='vertical', shrink=0.7)
# This makes the original colorbar look a bit out of place,
# so let's improve its position.
l, b, w, h = plt.gca().get_position().bounds
ll, bb, ww, hh = CB.ax.get_position().bounds
CB.ax.set_position([ll, b + 0.2*h, 1.2*ww, h*0.3])
#for i in range(len(label_ind)):
#    plt.text(label_x[i],label_y[i], label_ind[i])
texts = []
for x, y, text in zip(label_x, label_y, label_ind):
    texts.append(plt.text(x, y, text))
adjust_text(texts, precision=0,force_text=0.01, size=7,autoalign='y',only_move={'points':'x', 'text':'x'},color='r')
#adjust_text(texts, force_text=0.05, arrowprops=dict(arrowstyle="-|>",
#                                                    color='r', alpha=0.5))
#adjust_text(texts,prefer_move='y', precision=0, arrowprops=dict(arrowstyle="->", color='r', lw=0.5))
#CBI.ax.set_position([ll, b + 0.1*h, ww, h*0.9])
#TODO 
# 2. Get slices (see below, needs to be integrated with the above)
# 3. Add background
# 4. Include RLPs just beyond the field of view and fix code
#plt.show()
# 2d image plot with profiles
#plt.figure(2,figsize=(10,10)) # default seem to be centimeters
#h, w = data.shape
#gs = gridspec.GridSpec(2, 2,width_ratios=[w,w*.2], height_ratios=[h,h*.2])
#ax = [plt.subplot(gs[0]),plt.subplot(gs[1]),plt.subplot(gs[2])]
#bounds = [x.min(),x.max(),y.min(),y.max()]
#ax[0].imshow(data, cmap='gray', extent = bounds, origin='lower')
#ax[1].plot(data[:,w/2],Y[:,w/2],'.',data[:,w/2],Y[:,w/2])
#ax[1].axis([data[:,w/2].max(), data[:,w/2].min(), Y.min(), Y.max()])
#plt.plot(XI[h/2,:],ZI[h/2,:],'.',XI[h/2,:],ZI[h/2,:])
#plt.show()
plt.figure(2,figsize=(8.5,8.5)) # default seem to be centimeters
ZI=np.transpose(Z)
# 02 11 12
# for a 3x3x3
#hc = [0,0,0]; kc = [0,3,3]; lc = [6,3,6]
# for a 1x6x6
#hc = [0,0,0]; kc = [0,6,6]; lc = [12,6,12]
# for a 6x2x2
#hc = [0,0,0]; kc = [0,2,2]; lc = [4,2,4]
l0 = [];l1=[]
k = len(hc)
for i in range (k):
    Ghkl = hc[i]*astar+kc[i]*bstar+lc[i]*cstar
    Gperp = abs(np.dot(Ghkl,astarhat))
    Gmag2 = np.dot(Ghkl,Ghkl)
    Gpara = np.sqrt(abs(Gmag2-Gperp*Gperp))
#    print Gpara
    ind = int((Gpara/delta)-1)
    l0.append(ind)
    l1.append(ind+3)
scan0 = [] 
#i0 = int(1.32/delta)-1;i1 = int(1.56/delta)-1;i2 = int(1.89/delta)-1
if (store == 'true'):
    f = open(sfile,'w')
llmax = []
for j in range (k):
    scan0 = []
    for i in range (len(XI[j])):
        scan0.append(ZI[i,l0[j]:l1[j]].sum())
    llmax.append(np.nanmax(scan0))
    plt.plot(XI[j],scan0)
    if (store == 'true'):
        line = '#'+str(hc[j])+' '+str(kc[j])+' '+str(lc[j])+' \n'
        f.write(line)
        line = '# \n'
        f.write(line)
        k = 0
        for i in XI[j]:
#        print line
            line="%10.5f%12.1f" % (i,scan0[k])
            line = line + "\n"
            f.write(line)
            k += 1
#        print 'Not done yet but we need to store the plots for later use'
if (store == 'true'):
    f.close
# Now for experimental data
#lldata="/home/winokur/Documents/gixrd/NaT3/NaT3_ll_expt.dat"
#lldata='none'
#print maxcu
if (lldata != 'none'):
    xs=0.
    ys= 750.
    ys= 0.
#    maxcu = np.amax(llmax[1])
    maxcu = llmax[1]
#    print ' maxcu', maxcu
    ict,xd,yd0,yd1,yd2 = read_data_file(lldata)
    ym = np.array([1.,1.,1.,1.])
#    ys = np.array([15000.,10000.,5000.,0.])
# Now to plot the experimental layer lines
    xxd = np.empty((ict), dtype=object) # Define an empty array
    yyd = np.empty((ict), dtype=object) # Define an empty array

    for i in range(ict):
        xxd[i] = xd[i]*1.05+xs
    maxda = np.amax(yd2)
    ym[0] = maxcu/(maxda-ys)/2.8
    for i in range(ict):
        yyd[i] = (yd2[i]-ys)*ym[0]
        if (yyd[i]<0.0):
            yyd[i]=0.
    plt.plot(xxd,yyd, label = '$k=1$'+'$,\ell=1$',linewidth=1.0,linestyle='-.',color='black')
    for i in range(ict):
        yyd[i] = (yd1[i]-ys)*ym[0]
        if (yyd[i]<0.0):
            yyd[i]=0.
    plt.plot(xxd,yyd, label = '$k=0$'+'$,\ell=2$',linewidth=1.0,linestyle='--',color='black')
#    for i in range(ict):
#        yyd[i] = yd0[i]*ym[0]
#s    plt.plot(xxd,yyd,label = '$k=$'+'$,\ell=$',linewidth=1.0,linestyle='solid',color='black')

    plt.legend(loc='upper right',fontsize=12)  # do this after the labels and plots
    plt.xlabel(r'a* ($\rm{\AA}^{-1}$)',fontsize=14)
    plt.ylabel('Intensity (arb. units)',fontsize=14)


#if (k > 0):
plt.show()