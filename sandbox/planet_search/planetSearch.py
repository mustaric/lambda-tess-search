#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  9 16:47:22 2019

@author: smullally
"""

from astropy.timeseries import BoxLeastSquares
import numpy as np
from astropy.convolution import convolve, Box1DKernel


def clean_timeseries(time, flux, qflags, det_window, noise_window, n_sigma):
    
    qbad = qflags != 0
    bad = idNoisyData(flux, noise_window, Nsigma=n_sigma)
    #qbad and bad should be of the same length
    
    flagged = bad | qbad  #Indicate bad data
    med_det = median_detrend(flux, det_window, flagged)
    good_time = time[~flagged]
    
    return good_time, med_det

def median_detrend(flux,window, gap):
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
        
    mdFlux = median_detrend(flux, win, is_bad)

    for i in np.arange(1,4):

        if np.all(is_bad):
            continue

        sd = np.std(mdFlux[~is_bad])
        is_bad |= np.abs(mdFlux) > Nsigma * sd

    return is_bad

def is_odd(num):
   return num % 2 != 0


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
    """
    
    model = BoxLeastSquares(time,flux)
    fmodel =  model.model(time,results[0],results[3],results[1])
    noise = np.std(flux-fmodel)
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
    print(maxTces)
    keepLooking = True
    counter = 0
    results = []
    stats = []
    bls_durs_day=np.array(bls_durs_hrs)/24
    
    t=time.copy()
    f=flux.copy()
    
    #print(t[0:10])
    #print(f[0:10])
    
    while keepLooking:
        
        bls_results = findBlsSignal(t, f, bls_durs_day, minP=minP, maxP=maxP)
        #print(bls_results)
        bls_results[4] = simpleSnr(time,flux,bls_results)
        #simple ssnr because the BLS depth snr is acting strangely
        
        results.append(bls_results)
        bls = BoxLeastSquares(t,f)
        bls_stats = bls.compute_stats(bls_results[0], bls_results[3],bls_results[1])
        stats.append(bls_stats)
        #signal_snr = bls_stats['depth'][0]/bls_stats['depth'
        transit_mask = bls.transit_mask(t, bls_results[0],\
                                        bls_results[3]*1.1, bls_results[1])
        
        t=t[~transit_mask]
        f=f[~transit_mask]
        
        
        if (len(t)/len(time) > fracRemain) & \
               (bls_results[4] >= minSnr) & \
               (counter <= maxTces) :
                
            counter=counter + 1
            keepLooking = True
            
        else:          
            keepLooking = False
            print(len(t)/len(time))
            print(counter)
 

    return np.array(results), np.array(stats)




