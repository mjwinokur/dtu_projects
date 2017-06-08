#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 11:50:16 2017

@author: winokur
"""

#import numpy as np
import matplotlib.pyplot as plt

def energy_time_plot(tsteps,lj_E,dipl_E,bd_E,ang_E,opls_E,imp_E,totalE):
#def energy_time_plot(tsteps,lj_E,dipl_E,total_E):
# evenly sampled time at 200ms intervals
#t = np.arange(0., 5., 0.2)
# red dashes, blue squares and green triangles
    plt.plot(tsteps, lj_E, 'r--', tsteps, dipl_E, 'bs', tsteps, bd_E, 'g^',tsteps,ang_E,'b-',
             tsteps,opls_E,'g-',tsteps,totalE,'r-')
    plt.show()
    return
#
'''
tsteps =[]; lj_E =[]; dipl_E = []; bd_E = []; ang_E = []; opls_E = []; imp_E =[]; totalE = []
tsteps = np.arange(0.,20000.,2000.)
for i in range (len(tsteps)):
    lj_E.append(tsteps[i]/10.)
    dipl_E.append((tsteps[i]/1000.)**2.)
    bd_E.append(tsteps[i]/5.)
    opls_E.append(0.)
    ang_E.append(0.)
    imp_E.append(0.)
    totalE.append(lj_E[i]+dipl_E[i]+bd_E[i]+opls_E[i])
energy_time_plot(tsteps,lj_E,dipl_E,bd_E,ang_E,opls_E,imp_E,totalE)
print 'done'
'''