#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 13:27:27 2017

@author: winokur
"""
import numpy as np
from scipy import fftpack
import matplotlib.pyplot as plt

a = 0.5
N = 500 #number of points
x = np.linspace(-5,5,N+1)
x2 = np.linspace(-10,10,2*N+1)
lorentz = (a/np.pi) * (1/(a**2 + x2**2)) #lorentzian function
delta = x*0.0
delta[N/2]=1.
fourier = fftpack.rfft(lorentz)

fourier = (2 * np.pi**2 / N)* abs(fourier[0:N/4]) #supposed to normalize
mu=0.;sig=0.25;ab=1.0
gauss = ab*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
inv_area = 1./np.sum(lorentz)
test = inv_area*np.convolve(lorentz,gauss, mode='full')
#test = inv_area*np.convolve(lorentz,gauss, mode='valid')
#test = np.convolve(lorentz,delta, mode='valid')
plt.plot(x2,lorentz)
plt.plot(x,gauss)
m=len(test)/2
sub=test[m-N/2:m+N/2+1]
plt.plot(x,sub)
plt.show()
