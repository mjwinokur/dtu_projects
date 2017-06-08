#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 20:59:56 2017

@author: winokur
"""
import numpy as np
from mw_library import read_txyz
from mw_library import axisangle_to_q
#from mw_library import q_to_axisangle
from mw_library import qv_mult
import scipy.optimize
import functools

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

txyzname = "NaT2mw.txyz" # Base monomer in tinker format and includes the unit cell properties on line 2
l2_num,aax,aay,aaz,atype,astring,abt,uc,header = read_txyz(txyzname)
x = np.append(np.asarray(aax[0:24],dtype=float),np.asarray(aax[48:72],dtype=float))
y = np.append(np.asarray(aay[0:24],dtype=float),np.asarray(aay[48:72],dtype=float))
z = np.append(np.asarray(aaz[0:24],dtype=float),np.asarray(aaz[48:72],dtype=float))
"""
x = []; y = []; z=[]
for i in range (24):
    x = np.append(x,float(aax[i]))
    y = np.append(y,float(aay[i]))
    z = np.append(z,float(aaz[i]))
"""
data = np.concatenate((x[:,np.newaxis],y[:,np.newaxis], z[:,np.newaxis]),axis=1)
# Calculate the mean of the points, i.e. the 'center' of the cloud
datamean = data.mean(axis=0)
temp_dat = data - datamean
# Do an SVD on the mean-centered data.
uu, dd, vv = np.linalg.svd(temp_dat)
#
# Now vv[0] contains the first principal component, i.e. the direction
# vector of the 'best fit' line in the least squares sense.
#
# Shift by the mean to get the line in the right place
linepts = np.array([datamean,datamean+4.*vv[0]])
# Find best-fit plane
fun = functools.partial(error, points=data)
params0 = [0, 0, 0]
res = scipy.optimize.minimize(fun, params0)
a = res.x[0]; b = res.x[1]; c = res.x[2]
point  = np.array([0.0, 0.0, c])
norm_out = np.cross([1,0,a], [0,1,b])  #Normal to the plane
d = -point.dot(norm_out)
#xx, yy = np.meshgrid([-1,20], [-1,20])
#z = (-norm_out[0] * xx - norm_out[1] * yy - d) * 1. /norm_out[2]
norm_outP = np.array([datamean,datamean+4.*norm_out])
norm_in=np.cross(norm_out, vv[0])
norm_inP = np.array([datamean,datamean+4.*norm_in])

dataR =np.array(temp_dat)
q1R1 = axisangle_to_q(vv[0],np.pi) # 180 deg flip around long axis
q1R2 = axisangle_to_q(norm_in ,0.50*np.pi)
q1R3 = axisangle_to_q(norm_out,0.50*np.pi)
for i in range(48):
    temp = dataR[i]
    temp = qv_mult(q1R1,(temp[0],temp[1],temp[2])) 
    temp = qv_mult(q1R2,(temp[0],temp[1],temp[2])) 
    temp = qv_mult(q1R3,(temp[0],temp[1],temp[2])) 
    dataR[i]=temp
dataR = dataR + datamean
#
# Verify that everything looks right.
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as m3d
ax = m3d.Axes3D(plt.figure())
ax.set_xlim3d(-1,20)
ax.set_ylim3d(-1,20)
ax.set_zlim3d(-1,20)
ax.scatter3D(*data.T)
ax.scatter3D(*dataR.T)
ax.plot3D(*linepts.T)
#ax.plot_surface(xx, yy, z, alpha=0.2, color=[0,1,0])
ax.plot3D(*norm_outP.T)
ax.plot3D(*norm_inP.T)
plt.show()