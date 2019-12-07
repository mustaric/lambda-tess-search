#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 16:11:32 2019

Adapted from DAVE to remove the clipboard information.
Expects only a raw and detrended light curve.

@author: smullally
"""

import matplotlib.pyplot as plt
import numpy as np

def plotDiagnosticLightcurves(time, rawLc, detrendLc, path="."):

    plt.clf()
    #plt.gcf().set_size_inches((8,10))
    plt.subplot(211)
    plt.plot(time, rawLc, 'k.')
    plt.ylabel("Raw Lightcurve")


    plt.subplot(212)
    plt.plot(time, detrendLc, 'k.')
    plt.xlabel('Time (BTJD)')
    plt.ylabel("Detrended Lightcurve")


def plotData(time_raw, raw, time_det, detrend, meta, mPanel=3, nPanel=3):
    """Plot raw and detrended lightcurves
    meta is a dictionary with tic sector, period, epoch, dur
    time_raw is same length as raw
    time_det is same length as det
    oPanel is offset in panels *2
    """

    fl = np.zeros(len(time_raw)) == 1
    
    markTransits = True
    try:
        per = meta['period']
        epc = meta['epoch']
        dur_days = meta['dur']
    except KeyError:
        markTransits = False

    plt.gcf()
    #fig.set_size_inches(11, 8.5)

    colour = ["#FFF8F8", "#F8FFF8", "#F4F4FF"]
    start = np.min(time_raw[~fl])
    deltaT = np.max(time_raw[~fl]) - start
    deltaT /= float(nPanel)

    rawRange = np.percentile(raw[~fl], [1,99])

    for i in range(nPanel):
        ax = plt.subplot(2*mPanel*nPanel, 1, 2*i+1, facecolor=colour[i])
        plt.plot(time_raw[~fl], raw[~fl], 'ko', ms=2, alpha=.8)
        plt.ylim(rawRange)
        plt.ylabel("Raw flux")

        if markTransits:
            plotTransitRegions(time_raw[~fl], per, epc, dur_days)

        plt.subplot(2*mPanel*nPanel, 1, 2*i+2, sharex=ax, facecolor=colour[i])
        #Plotting bad data cadences turned off
#        plt.plot(time[fl], 0*time[fl], 'mo', ms=8, mec="none")
        plt.plot(time_det, detrend, 'ko', ms=2, alpha=.8)
        plt.ylabel("Detrended")


        if markTransits:
            plotTransitRegions(time_det, per, epc, dur_days)

        plt.xlim(start + i*deltaT, start + (i+1)*deltaT)
        plt.xlabel("Time (BTJD)")

    


def plotTransitRegions(time, period_days, epoch, duration_days, **kwargs):
    tmin = np.min(time)
    tmax = np.max(time)

    n1 = int(np.floor( (tmin-epoch)/period_days) )
    n2 = int(np.ceil(  (tmax-epoch)/period_days) )
    color = kwargs.pop('color','#888888')
    alpha = kwargs.pop('alpha', 1)

    for n in range(n1, n2+1):
        t0 = epoch + n*period_days
        lwr = t0 - .5*duration_days
        upr = t0 + .5*duration_days
        plt.axvspan(lwr, upr, color=color, alpha=alpha)

def plotFolded(time_det, detrend, model, meta, doublePeriod = False, modelOn = True):
    #plt.figure()  #Does not belong here. Remove from this function!
    #fl = np.zeros(len(time_det)) == 0
    period = meta['period']
    epoch = meta['epoch']
    ticid = meta['id']  #In prepartion for the multi-search pipeline
    flux = detrend
    time = time_det

    if doublePeriod:
        period *= 2
    if epoch > time[0]:
        diff = epoch-time[0]
        epoch -= period*np.ceil(diff/period)
        print("Reducing epoch")

    #plt.cla()
    phi = np.fmod(time-epoch + .25*period, period)

    plt.plot(phi, 1e6*flux, 'ko', ms=4)
    plt.plot(period+phi, 1e6*flux, 'o', color='#888888', ms=5, mec="none")
    #plt.title("ID: %s"%ticid)
    if modelOn:
        x = phi
        y = 1e6*model
        idx = np.argsort(x)
#
        plt.plot(x[idx], y[idx], 'r-', lw=1)
        plt.plot(period+x[idx], y[idx], 'r-')
    else:
        (a,b)=plt.ylim()
        plt.plot([.25*period,.25*period],[a,b],'r--')

    plt.ylabel("Fractional Amplitude (ppm)")
    plt.xlabel("Phase (days)")



def summaryPlot1(time_raw, raw, time_det, detrend, model, ave_image, meta):
    """
    Plot the data using the above functions.
    meta contains period, epoch, id, sector, dur, snr, depth, Tmag, Tstar
    time array
    detrend array
    model array of transit model
    output is fill name to write the plot to.
    """
    ticid=meta['id']
    snr=meta['snr']
    per=meta['period']
    dur=meta['dur']
    depth=meta['depth']
    pn = meta['pn']
    imloc = meta['imloc']  # (col,row)
    radius = meta['radius']
    sector = meta['sector']
    cam = meta['cam']
    ccd= meta['ccd']
    
    plt.clf()
    #plt.subplot(3,2,(1,2))
    plotData(time_raw, raw, time_det, detrend, meta, nPanel=1)
    titlewords="TIC=%s pn:%i Sector:%i, Cam: %1u, CCD: %1u\nP=%.4f days, Dur=%.2f hr, Depth=%.1f ppm, SNR=%.2f " \
        % (ticid, pn, sector, cam, ccd, per, dur*24, depth * 1e6,snr)
    plt.annotate(titlewords,(.1,.92),xycoords="figure fraction", fontsize=14)
    #plt.suptitle("TIC %i Sector %i" %(meta['id'], meta['sector']))
    
    plt.subplot(3,2,(3,4))
    plotFolded(time_det, detrend, model, meta, modelOn=True)
    
    plt.xlim((0,per))


    plt.subplot(325)
    plotFolded(time_det, detrend, model, meta, modelOn=True)
    if dur*1.5 < per:
        plt.xlim(per*.25-dur*1.9,per*.25+dur*1.9)
    else:
        plt.xlim(per*.15, per*.35)
    
    plt.subplot(326)
    plt.imshow(ave_image, vmax=np.percentile(ave_image, 96), cmap='cividis',origin='lower')
    plt.plot(imloc[0], imloc[1],'o',mfc="None", mec='red',mew=2, alpha=.3,ms=11)
    #circ=plt.Circle(imloc, radius, fill=False, lw=2)
    plt.title('Radius: %.2f px' % radius)
    
