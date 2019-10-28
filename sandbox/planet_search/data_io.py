#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 14:44:55 2019

@author: smullally
"""
import boto3
from astropy.io import fits
import os
from astropy.table import Table
import numpy as np


def read_header(phead):
    ticid = phead['OBJECT'][4:]
    camera = phead['CAMERA']
    sector = phead['SECTOR']
    ccd = phead['CCD']
    
    return ticid, camera, sector, ccd

def read_lightcurve_tasoc(bucket, bucket_filename, local_filename):
    
    s3 = boto3.resource('s3')
    bucket = s3.Bucket(bucket)
    
    bucket.download_file(bucket_filename, local_filename, 
                         ExtraArgs={"RequestPayer": "requester"})
    data = fits.getdata(local_filename, 1)
    
    time = data['TIME']
    flux = data['FLUX_RAW']
    qflags = data['PIXEL_QUALITY']
    
    os.remove(local_filename)
    
    return time, flux, qflags


def read_lightcurve(filename):
    data = fits.getdata(filename, 1)
    head = fits.getheader(filename, 1)
    time = data['TIME']
    flux = data['SAP_FLUX']
    qflags = data['QUALITY']
    
    return time, flux, qflags, head

def read_lightcurve_lambda(bucket, bucket_filename, local_filename):
    """
    read in lambda generated light curves.
    """

    s3 = boto3.resource('s3')
    bucket = s3.Bucket(bucket)
    
    bucket.download_file(bucket_filename, local_filename, 
                         ExtraArgs={"RequestPayer": "requester"})
    
    time,flux, qflags, head = read_lightcurve(local_filename)
    
    os.remove(local_filename)
    
    return time, flux, qflags, head

def read_lightcurve_lambda_local(bucket, bucket_filename, local_filename):
    """
    read in lambda generated light curve from local directory
    For testing.
    In this case, only local_filename matters. It is your data.
    The code below should match the code above.
    """

    time,flux, qflags, head = read_lightcurve(local_filename)
    
    #os.remove(local_filename)
    
    return time, flux, qflags, head

def write_results(filename, ticid, results, stats, **kwargs):
    
    header ='bls_period_day,bls_t0_btjd, bls_depth_frac,\
            bls_duration_day, bls_snr, n_trans_good, tic.pn'
    info_array=np.zeros([len(stats),7])
    #info_array[:,6] = ticid
    
    for i,r in enumerate(results):
        astat = stats[i]
        ntrans = getNumberTransits(astat, threshold=1e-12)
        info = np.concatenate((r,[ntrans]))

        info_array[i,:-1] = info.transpose()
        info_array[i,6] = ticid + (i+1)/100.0
        
    hdrStr = make_header(**kwargs) + header
    np.savetxt(filename, info_array, newline = '\n', delimiter = ',', header = hdrStr)
    
def make_header(**kwargs):
    output = []
    for k in kwargs.keys():
        output.append("%s: %s" %(k, kwargs[k]))
        
    return "\n".join(output) + "\n"

def getNumberTransits(astat, threshold = 1e-12):
        good = astat['per_transit_log_likelihood'] >= threshold
        ntransits = len(good[good])
        
        return ntransits

def write_timeseries(filename, time, flux, head):
    
    #Need to fix this to pass some of he header information along
    t = Table([time, flux], names=('TIME','DETREND_FLUX'))
    t.write(filename, format='fits', overwrite=True)