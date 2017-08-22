#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 10:25:37 2017

@author: winokur
Function to determine the angle beween two planes formed by four atoms with the 2nd position as the central atom
"""
import numpy as np
#import math
#import os
#from mayavi.mlab import *
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
#
def improper_calc_ang(p1, p2, p3, p4):
    va1 = np.array([p1[0]-p2[0],p1[1]-p2[1],p1[2]-p2[2]])
    va2 = np.array([p3[0]-p2[0],p3[1]-p2[1],p3[2]-p2[2]])
    vperpa=np.cross(va1,va2)
    vb1 = va2
    vb2 = np.array([p4[0]-p2[0],p4[1]-p2[1],p4[2]-p2[2]])
    vperpb=np.cross(vb1,vb2)
    return angle_between(vperpa,vperpb)
#
"""
p1 = [7.2787, 2.7807, 0.9833]
p2 = [8.322, 3.9043, 0.2102]
p3 = [9.645, 2.8165, 0.2588]
p4 = [9.2798, 1.623, 0.8561]
angle = improper_calc_ang(p2,p1,p3,p4)
print angle
"""