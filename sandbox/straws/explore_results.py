#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 15:29:16 2019

@author: smullally
"""

import boto3
import numpy as np
from astropy.timeseries import BoxLeastSquares
import matplotlib.pyplot as plt


 #%%   
def plot_bls_folds(time,flux,results,ticid, zoom=True):
    
    num=4  #max number of possible planet candidates
    plt.figure(figsize=(8,12))
    s=np.std(flux)
    
    for i,r in enumerate(results):
        plt.subplot(num,1,i+1)
        phase = (time - r[1] + 0.5*r[0]) % r[0] - 0.5*r[0]
        
        model = BoxLeastSquares(time,flux)
        fmodel = model.model(time,results[i,0],results[i,3],results[i,1])
        order=np.argsort(phase)
        plt.plot(phase[order]*24,fmodel[order],'g.-', ms=3)
        plt.plot(phase*24,flux,'k.',label="TIC %u P=%6.2f d" % (ticid, results[i,0]))
        #plt.title(str(ticid) + " Per: " + str(results[i,0]) )
        if zoom:
            plt.xlim(-5.7*results[i,3]*24,5.7*results[i,3]*24)
            plt.ylim(-2.9*r[2], 4*s)
        
        plt.legend(fontsize=8)




#%%

search_bucket = "tesssearchresults"
lc_bucket = "straw-lightcurves"

ticid = 126944572
camera = 1
ccd = 1
sector = 1


search_path = "tic%012u/" % ticid
rootname = "tic%012u_s%04u-%1u-%1u" %  (ticid, sector, camera, ccd)
bls_name = "%s_plsearch.csv" % rootname
det_name = "%s_detrend.fits" % rootname
lc_name = "%s_stlc.fits" % rootname

print(bls_name,det_name,lc_name)

client = boto3.client("s3")




        