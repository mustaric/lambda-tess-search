#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 30 13:25:47 2020

@author: smullally
"""


import numpy as np
import planetSearch as ps
import gen_lightcurve as lc
import matplotlib.pyplot as plt
from random import gauss

#Use one star as an example.
#Then produce Gaussian noise ight curve at the same time stamps 

ticid = 350814695	#(near to kepler10 but not kepler 10, in S15)
sector = 15
time, flux, flags = lc.eleanor_corr(ticid, sector)
#%%
det_window = 23
noise_window = 15
n_sigma = 4  #noise reject sigma

good_time, meddet_flux = ps.clean_timeseries(time, flux, flags, det_window, \
                                          noise_window, n_sigma, sector)

    #%%
plt.figure()
plt.subplot(312)
plt.plot()
plt.plot(good_time, meddet_flux,'.')
plt.subplot(311)
plt.plot(time,flux,'.')

#%%
lcstd = np.std(meddet_flux)

rand_flux = np.array([gauss(0.0, lcstd) for i in range(len(meddet_flux))])

plt.subplot(313)
plt.plot(good_time,rand_flux,'.')

#%%
results, stats = ps.identifyTces(good_time, rand_flux, bls_durs_hrs=[1,2,4,8,12],
                                 minSnr=.5, fracRemain=0.5, 
                                 maxTces=30, minP=None, maxP=7.5)
print(results)
#results[4] contains an snr
#results = period, t0, depth, duration, snr
#%%
Prange = [.4, 7] #period range in days
minSnr = 3
fracRemain = 0.45
bls_durs = [4]
Nlcs = 2000
maxTces = 1
all_results = np.zeros([Nlcs,5])
lcstd = 0.00316  #facot of 10 larger than usual


for i in range(Nlcs):
    
    rand_flux = np.array([gauss(0.0, lcstd) for i in range(len(meddet_flux))])
    
    results, stats = ps.identifyTces(good_time, rand_flux, bls_durs_hrs=bls_durs,
                                 minSnr=minSnr, fracRemain=fracRemain, 
                                 maxTces=maxTces, minP=Prange[0], maxP=Prange[1])
    
    all_results[i,:] = results[0,:]
    

#results1 - lowsnr
#results2 = 10x snr

