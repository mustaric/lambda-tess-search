#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  9 16:47:22 2019

@author: smullally
"""

from astropy.timeseries import BoxLeastSquares
import numpy as np
from astropy.convolution import convolve, Box1DKernel

#import matplotlib.pyplot as plt


def clean_timeseries(time, flux, qflags, det_window, noise_window, n_sigma):
    
    qbad = qflags != 0
    bad = idNoisyData(flux, noise_window, Nsigma=n_sigma)
    #qbad and bad should be of the same length
    
    flagged = bad | qbad  #Indicate bad data
    med_det = median_detrend(flux[~flagged], det_window)
    det_time = time[~flagged]
    
    #Look for 3 bad sections on length of around 2 days (window = 90)
    std_bad, med_std = running_std_gap(med_det, 70, N=3, ntimes=3)
    #print(len(std_bad[std_bad]))
    
    
    good_time = det_time[~std_bad]
    good_flux = med_det[~std_bad]
    
    #plt.figure()
    #plt.plot(time,flux/np.nanmedian(flux) - 1, 'ko')
    #plt.plot(good_time, good_flux, 'r.')
    
    return good_time, good_flux

def median_detrend(flux, window):
    """
    Fergal's code to median detrend. 
    """
    size = len(flux)
    nPoints = window
    
    filtered = np.zeros(size)
    for i in range(size):
        #This two step ensures that lwr and upr lie in the range [0,size)
        lwr = max(i-nPoints, 0)
        upr = min(lwr + 2*nPoints, size)
        lwr = upr- 2*nPoints

        sub = flux[lwr:upr]

        offset = np.median(sub)
        try:
            filtered[i] = flux[i]/offset - 1
        except ZeroDivisionError:
                filtered[i] = 0

    return filtered

def median_subtract(flux, window):
    """
    Fergal's code to median detrend. 
    """
    size = len(flux)
    nPoints = window
    
    filtered = np.zeros(size)
    for i in range(size):
        #This two step ensures that lwr and upr lie in the range [0,size)
        lwr = max(i-nPoints, 0)
        upr = min(lwr + 2*nPoints, size)
        lwr = upr- 2*nPoints

        sub = flux[lwr:upr]

        offset = np.median(sub)
        try:
            filtered[i] = flux[i] - offset
        except ZeroDivisionError:
                filtered[i] = 0

    return filtered

def conv_detrend(flux, window, gap):
    """
    return median detrended array
    centered on zero at the end.
    if gap is true, then don't use that point
    """

    #filt = medfilt(flux,window) + 1e-8
    box_kernel = Box1DKernel(window, mode="linear_interp")
    filt = convolve(flux[~gap], box_kernel, boundary="extend")
    median_det = flux[~gap]/filt - 1
    
    return median_det

def idNoisyData(flux, window, Nsigma=4):
    """
    Determine sections of the data that are very noisy compared to the rest.
    Look for N sigma away from the std of the data.
    Be careful, try not to cut out planets.
    I recommend using a window that is smaller than used for planet finding.
    
    """
    win = int(window)
    if ~is_odd(win):
        win=win+1

    is_bad = np.zeros(len(flux)) == 1
        
    mdFlux = median_detrend(flux[~is_bad], win)

    for i in np.arange(1,4):

        if np.all(is_bad):
            continue

        sd = np.std(mdFlux)
        is_bad |= np.abs(mdFlux) > Nsigma * sd

    return is_bad

def is_odd(num):
   return num % 2 != 0


def running_std_gap(flux, window, N=3, ntimes=3):
    """
    for specified window, determine data chunks that are parts of sections
    of the data that have std ntimes larger than the overall std. only pulls
    out N sections.
    
    Returns isbad array of 1 and 0 where 1 means bad data and 0 means clean
    Be sure to set a wide enough window so you don't throw away planets.
    Probably N*duration(in points) of longest transit expected.
    """
    gap = np.zeros(len(flux))
    
    std_array = np.zeros(len(flux))
    
    
    for i in range(window, len(flux), 1):
        
        f = flux[i-window:i]
        std_array[i] = np.nanstd(f)
    
    #plt.figure()
    #plt.plot(std_array,'.')
    
    med_std = np.median(std_array)
    #print("mean std: %f" % med_std)
    
    argsort = np.argsort(std_array)
    
    for i in argsort[-1*N:]:
        if std_array[i] > med_std * N:
            gap[i-window:i] = 1
            
    isbad = gap == 1
    
    return isbad, med_std
        

def calcBls(flux,time, bls_durs, minP=None, maxP=None, min_trans=3):
    """
    Take a bls and return the spectrum.
    """
    
    bls = BoxLeastSquares(time, flux)
    period_grid = bls.autoperiod(bls_durs,minimum_period=minP, \
                   maximum_period=maxP, minimum_n_transit=min_trans, \
                   frequency_factor=0.8)
    
    bls_power = bls.power(period_grid, bls_durs, oversample=20)
    
    return bls_power

def findBlsSignal(time, flux, bls_durations, minP=None, maxP=None, min_trans=3):
    
    bls_power = calcBls(flux, time, bls_durations, minP=minP, maxP=maxP, min_trans=min_trans)
    
    index = np.argmax(bls_power.power)
    bls_period = bls_power.period[index]
    bls_t0 = bls_power.transit_time[index]
    bls_depth = bls_power.depth[index]
    bls_duration = bls_power.duration[index]
    bls_snr = bls_power.depth_snr[index]
    
    return np.array([bls_period, bls_t0, bls_depth, bls_duration, bls_snr])
 
def simpleSnr(time,flux,results):
    """
    calculate a simple snr on the planet based on the depth and scatter
    after you remove the planet model from the data.
    """
    
    model = BoxLeastSquares(time,flux)
    fmodel =  model.model(time,results[0],results[3],results[1])
    flatten = median_subtract(flux-fmodel, 12)
    
    noise = np.std(flatten)
    snr = results[2]/noise
    
    return snr

   
def identifyTces(time, flux, bls_durs_hrs=[1,2,4,8,12], minSnr=3, fracRemain=0.5, \
                 maxTces=10, minP=None, maxP=None):
    """
    Find highest point in the bls.
    remove that signal, median detrend again
    Find the next signal.
    Stop when less than half the original data set remains.
    Or, when depth of signal is less than snr*running_std 
    
    returns period, t0, depth, duration, snr for each signal found.
    """
    
    keepLooking = True
    counter = 0
    results = []
    stats = []
    bls_durs_day=np.array(bls_durs_hrs)/24
    
    t=time.copy()
    f=flux.copy()
    

    while keepLooking:
        
        bls_results = findBlsSignal(t, f, bls_durs_day, minP=minP, maxP=maxP)
        #print(bls_results)
        #simple ssnr because the BLS depth snr is acting strangely
        bls_results[4] = simpleSnr(t, f, bls_results)
        
        
        results.append(bls_results)
        bls = BoxLeastSquares(t,f)
        bls_stats = bls.compute_stats(bls_results[0], bls_results[3],bls_results[1])
        stats.append(bls_stats)
        #signal_snr = bls_stats['depth'][0]/bls_stats['depth'
        transit_mask = bls.transit_mask(t, bls_results[0],\
                                        bls_results[3]*1.1, bls_results[1])
        #plt.figure()
        #plt.plot(t,f,'ko',ms=3)
        
        t=t[~transit_mask]
        f=f[~transit_mask]
        
        #plt.plot(t,f,'r.')
        #Conditions to keep looking
        if (len(t)/len(time) > fracRemain) & \
               (bls_results[4] >= minSnr) & \
               (counter <= maxTces) :
                
            counter=counter + 1
            keepLooking = True
            
        else:          
            keepLooking = False
 

    return np.array(results), np.array(stats)




