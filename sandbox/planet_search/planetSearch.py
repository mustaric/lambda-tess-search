#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  9 16:47:22 2019

@author: smullally
"""

from astropy.timeseries import BoxLeastSquares
import numpy as np
from astropy.convolution import convolve, Box1DKernel
from plateau import plateau

#import matplotlib.pyplot as plt


def clean_timeseries(time, flux, qflags, det_window, noise_window, n_sigma, sector):
    
    qbad = qflags != 0
    tgaps = loadGapInfoBySector(time,sector)
    bad = idNoisyData(flux, noise_window, Nsigma=n_sigma)
 
    flagged = bad | qbad | tgaps  #Indicate bad data
    med_det = median_detrend(flux[~flagged], det_window)
    det_time = time[~flagged]
    
    #Look for 3 bad sections on length of around 2 days (window = 90)
    std_bad, med_std = running_std_gap(med_det, 70, N=3, nSigTimes=4.5)
    #print(len(std_bad[std_bad]))
    
    
    good_time = det_time[~std_bad]
    good_flux = med_det[~std_bad]
    
    #Final Pass for single point outliers that are not periodic
    spo_idx = findOutliers(good_time, good_flux, gap=None,
                 threshold_sigma=2.75,
                 precision_days=0.0209,
                 maxClusterLen = 2
                 )   
    #print("spo_idx")
    #print(spo_idx)
    #spogaps = np.zeros(len(spo_idx)) == 1
    #spogaps[spo_idx] = True
    #plt.figure()
    #plt.plot(good_time[~spo_idx], good_flux[~spo_idx])
    
    return good_time[~spo_idx], good_flux[~spo_idx]

def loadGapInfoBySector(time, sector):
    """Loads a list of bad cadence indices.

    TESS produces quality flags, but  does not populate the FFIs with them.
    Instead we have to look the up in the data release notes.

    Based on the Data release notes, but modified by hand based on
    inspection of Wasp 126

    Inputs
    ---------
    time
        (1d np array) Array of TJDs for the data. See `extractlc.loadSingleSector`.
    sector
        (int)

    Returns
    -----------
    1d boolean np array of length time
    """
    gaps = np.zeros_like(time, dtype=bool)

    if sector == 1:
        # See page 2 of
        #https://archive.stsci.edu/missions/tess/doc/tess_drn/tess_sector_01_drn01_v01.pdf
        gaps |= (time <= 1325.61)
        gaps |= (1338.52153 <= time) & (time <= 1339.65310)  #Inter orbit gap
        gaps |= (1346.95 <= time) & (time <= 1349.75)  #See DRN 1 p3
        gaps |= (time >= 1352.92)  #End of sector usually bad.
    elif sector == 2:
        gaps |= (1367.15347 <= time) & (time <= 1368.59406)  #
    elif sector == 3:
        gaps |= (1381.1 <= time) & (time <= 1385.89663)  #Pre "science start"
        gaps |= (1394.47997 <= time) & (time <= 1395.80497)  #apears to be bad??
        gaps |= (1395.47997 <= time) & (time <= 1396.60497)  #Inter orbit gap
        gaps |= (1406.2 <= time) & (time <= 1409.38829)  #Post science
    elif sector == 4:
        #The bad guide star data may still be usable. Need to check
        gaps |= (1410.89974 <= time) & (time <= 1413.26468)  #Bad Guide star
#        gaps |= (1418.53691 <= time) & (time <= 1421.86)  #Instr. Anom.
#        gaps |= (1422.95 <= time) & (time <= 1424.54897)  #Inter orbit gap
        gaps |= (1418.53691 <= time) & (time <= 1424.54897)
        gaps |= (1436.0 <= time) & (time <= 1439.8)  #Not sure what this is
    elif sector == 5:
        gaps |= (1450.01 <= time) & (time <=  1451.81)  #Inter orbit gap
        #I don't think this means the data is generally bad
        gaps |= (1463.55 <= time) & (time <= 1464.40056)  #Camera 1 guiding
    elif sector == 6:
        gaps |= (1477.0 <= time) & (time <= 1478.41)  #inter orbit gap
        gaps |= (1463.6 <= time) & (time <= 1468.26998)  #before beginning of 6
    elif sector == 7:
        gaps |= (1517 <= time) & (time <= 1491.62) # orbit range
        gaps |= (1502.5 <= time) & (time <= 1505.01) #inter sector gap
    elif sector == 16:
        gaps |= (1738.65 <= time) & (time <= 1763.31) # orbit range
        gaps |= (1750.35 <= time) & (time <= 1751.652) #inter sector gap
        

#        gaps |= (<= time) & (time <= )  #

    else:
        raise ValueError("No gap info available for sector %i" %(sector))

    return gaps
    
    

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


def running_std_gap(flux, window, N=3, nSigTimes=3.3):
    """
    for specified window, determine data chunks that are parts of sections
    of the data that have std nSigTimes larger than the overall std. only pulls
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
        if std_array[i] > med_std * nSigTimes:
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


def findOutliers(time, flux, gap=None,
                 threshold_sigma=4,
                 precision_days=0.0205,
                 maxClusterLen = 2
                 ):
    """
    Identify single point outliers.

    Preserves consecutive outliers, and those that are evenly spaced in
    time. This protects short duration transits.

    Inputs:
    ------------
    time, flux
        (1d numpy array) Input data. Flux should have mean (or median) value
        of zero. Units of time are assumed to be days.

    Optional Inputs
    ------------------
    precision_days
        (float) Points that are evenly spaced to within this precision
        are considered periodic, and not marked as outliers. Setting
        this to zero turns off the search of periodicity.

    threshold_sigma
        (float) Points more than this many sigma from zero are considered
        potential outliers.
    maxClusterLen
        (int) Outliers are not marked if they are part of
        a contiguous cluster at least this long.

    Returns:
    ------------
    An array of indices indicating which points are single point
    outliers. The length of this array is equal to the number of outlier
    points flagged. The length of this array is NOT equal to the length
    of the flux array.

    Notes
    ----------
    `precision_days` should be set to a value comparable to the cadence time.
    For Kepler long cadence, this is 29.52 minutes = 0.0205 days.

    If `time` is not in units of days, set the value of precision in the
    same units.
    """

    assert not np.all(gap), "Can't find outliers if all data is gapped"

    if gap is None:
        gap = np.zeros_like(flux, dtype=bool)
    indices = np.zeros_like(gap)

    #Remove as much signal as possible from the data
#    fluxDetrended = medianDetrend(flux, 3)
    fluxDetrended = np.diff(flux)
    fluxDetrended = np.append(fluxDetrended, [0])  #Keep the length the same
    assert len(fluxDetrended) == len(flux)

    #Find outliers as too far away from the mean.
    rms = robustStd(fluxDetrended[~gap])
    threshold_counts = threshold_sigma * rms / np.sqrt(2)
    rawOutliers = plateau(np.fabs(fluxDetrended), threshold_counts)


    if len(rawOutliers) == 0:
        return indices

    # Remove clusters of 2 or more consectutive outliers
    #singleOutlierIndices = np.sort(outliers[(outliers[:,1] - outliers[:,0] <= 2)][:,0])
#    debug()
    span = rawOutliers[:,1] - rawOutliers[:,0]
    outliers = rawOutliers[span < maxClusterLen+2]
    for p1, p2 in outliers:
        indices[p1+1 :p2] = True


    #Check for periodicities in outliers
    if precision_days > 0:
        notOutliers = findPeriodicOutliers(time, indices, precision_days)
        indices[notOutliers] = False

    return indices

#from pdb import set_trace as debug

def findPeriodicOutliers(time_days, singleOutlierIndex, precision_days):
    """Identify groups of outliers that are equally spaced in time

    Inputs
    ---------
    time_days
        (1d numpy array) Times of data points, in days

    singleOutlierIndex
        (array of ints) Array elements of time that have
        been flagged as outliers. ``len(singleOutlierIndices == len(time)``

    precision_days
        (float) How close to perfectly evenly spaced do points need to
        be to be marked as periodic?

    Returns
    ----------
    notOutliers
        An array of elements of `singleOutlierIndices` that are periodic.
        set(notOutliers) is a wholly owned subset of set(singleOutlierIndices)
    """

    assert len(time_days) == len(singleOutlierIndex)

    #Convert to list of indices
    singleOutliers = np.where(singleOutlierIndex)[0]

    outlierTimes = time_days[singleOutliers]
    diffs = [outlierTimes[i+1] - outlierTimes[i] for i in range(0, len(outlierTimes)-1)]
    diffs = [round(d, 5) for d in diffs]

    if len(singleOutliers) >= 4:
        if len(set(diffs)) == len(diffs):
            possibleTimes = np.array([])
        else:
            period = max(set(diffs), key = diffs.count) # period = most common difference
            epoch = outlierTimes[ diffs.index(period) ]
            possibleTimes = np.arange(epoch, outlierTimes[-1] + 0.5*period, period)

        notOutliers = []
        for i in range(len(outlierTimes)):
            if np.any((abs(possibleTimes - outlierTimes[i]) < precision_days)):
                notOutliers.append(singleOutliers[i])
    

    elif len(singleOutliers) == 3:
        #If we only have three outliers, and they are equally spaced
        #then they are periodic
        if abs(diffs[0] - diffs[1]) < precision_days:
            notOutliers.extend(singleOutliers)
    #debug()
    assert set(notOutliers).issubset(singleOutliers)
    return notOutliers


def robustStd(y):
    assert len(y) > 0
    mad = y - np.median(y)
    mad = np.median(np.fabs(mad))

    #See https://en.wikipedia.org/wiki/Median_absolute_deviation
    #Valid for gaussians
    std = mad * 1.4826
    return std