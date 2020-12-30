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
#import lightkurve as lk
from astropy.io import fits
import io
import os
import pandas as p
import onepageReport as report

def generate_plots_s3(ticid,sector,cam,ccd, \
                   bls_bucket="tesssearchresults", \
                   detrend_bucket="tesssearchresults", \
                   ffilc_bucket="straw-lightcurves", \
                   outpath="/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/s0001-kdwarf/plots/"):
    
    """
    Given a TICID, sector, camera, ccd (all integers) and the bucket locations.
    generate a one page plot for each signal found by the bls using 
    the information in those files.
    bls_bucket contains the search results csv file.
    detrend_bucket contains the detrend fits file.
    ffi_bucket contains the raw light curve fits files.
    returns hdus and a dictionary of the csv.
    """
    
    rootname = "tic%012u_s%04u-%1u-%1u" %  (ticid, sector, cam, ccd)
    path = "tic%012u/" % ticid
    bls_name = "%s_plsearch.csv" % rootname
    det_name = "%s_detrend.fits" % rootname
    lc_name = "%s_stlc.fits" % rootname
    
    if detrend_bucket[0]=="/":
        det_hdu = fits.open(detrend_bucket+path+det_name)
        bls = np.loadtxt(bls_bucket+path+bls_name,delimiter=',')
    else:
        print(detrend_bucket,path,det_name)
        det_hdu = loadFitsFromUri(detrend_bucket, path, det_name)
        bls = loadCsvFromUri(bls_bucket, path, bls_name)
    if ffilc_bucket[0]=="/":
        raw_hdu = fits.open(ffilc_bucket+path+lc_name)
    else:
        raw_hdu = loadFitsFromUri(ffilc_bucket,path, lc_name)
    
    
    #Get the number of signals found by the bls.
    if len(bls.shape) == 1:
        N=1
        bls = bls.reshape((1,7))
        #print(bls.shape)
    else:
       N=bls.shape[0]

    
    time_raw = raw_hdu[1].data['TIME']
    raw = raw_hdu[1].data['SAP_FLUX']
    bkg = raw_hdu[1].data['SAP_BKG']
    time_det = det_hdu[1].data['TIME']
    detrend = det_hdu[1].data['DETREND_FLUX']
    """
    plt.figure()
    ax = plt.subplot(211)
    plt.plot(time_raw,raw,'.')
    plt.subplot(212, sharex=ax)
    plt.plot(time_raw,bkg, 'r.')
    """
    head = raw_hdu[1].header
    try:
        ave_im = raw_hdu[2].data
    except:
        ave_im = np.zeros((10,10))    
    meta={}
    meta['sector'] = sector
    meta['cam'] = cam
    meta['ccd'] = ccd
    try:
        meta['imloc'] = (head['CUBECOL'], head['CUBEROW'])
    except:
        meta['imloc'] = (head['APCEN_Y'], head['APCEN_X'])
        
    meta['radius'] = head['AP_RAD']
    
    
    for i in range(N):
        #print(i)

        meta['period'] = bls[i,0]
        meta['dur'] = bls[i,3]
        meta['epoch'] = bls[i,1]
        meta['snr'] = bls[i,4]
        meta['depth'] = bls[i,2]
        meta['ntrans'] = bls[i,5]
        meta['id'] = ticid
        meta['pn'] = i+1

       #print(meta)

               
        bls_object = BoxLeastSquares(time_det, detrend) 
        model = bls_object.model(time_det, meta['period'], \
                                      meta['dur'], meta['epoch'])
        out_name = "%s-%02i_plot.png" % (rootname, meta['pn'])
        output = outpath+out_name
        plt.figure(figsize=(10,12))
        report.summaryPlot1(time_raw, raw, time_det, detrend, model, ave_im, meta)
        
        plt.savefig(output)
        print(output)
    
    return raw_hdu[1].data

    
def loadFitsFromUri(bucketName, strawPath, fn):
    #boto stuff goes here
    s3=boto3.resource('s3') 
    uri = os.path.join(strawPath, fn)
    #print(uri)
    obj = s3.Object(bucketName, uri)
    thebytes = obj.get()['Body'].read()
    return fits.open(io.BytesIO(thebytes))    

def loadCsvFromUri(bucketName, strawPath, fn):
    #boto stuff goes here
    s3=boto3.resource('s3') 
    uri = os.path.join(strawPath, fn)
    #print(uri)
    obj = s3.Object(bucketName, uri)
    thebytes = obj.get()['Body'].read()
    return np.loadtxt(io.BytesIO(thebytes),delimiter=',') 


#%%
#December 9, 2019
#For S0001 run of Tmag 11-14 stars. all 60,000
path = "/Users/smullally/TESS/lambdaSearch/s01StrawRun60/blsResults"

good_ones = []
filenames = []

for root,dirs,files in os.walk(path):
    for f in files:
        if f[-10:] == "search.csv":
            blsResults = np.loadtxt(root+"/"+f, comments='#',delimiter=',')
            if len(blsResults) == 7:
                blsResults = np.reshape(blsResults,(1,7))
            for result in blsResults:
                if ((result[4]>4) & (result[4]<=7) & (result[5]>3)):
                    #print(result)
                    good_ones.append(result)
                    filenames.append(f)


#%%
#December 9, 2019
#A continuation of S0001 run. Run through good ones and make plots.
import matplotlib
matplotlib.use('Agg')
plotpath = "/Users/smullally/TESS/lambdaSearch/s01StrawRun60/blsResults/plot/second_tier/"
for filename in filenames:
    ticid = int(filename[3:15])
    sector = int(filename[20:21])
    cam = int(filename[22:23])
    ccd = int(filename[24:25])
    if sector == 1:    
        try:
            data = generate_plots_s3(ticid, sector,cam,ccd, outpath=plotpath)
        except NoSuchKey:
            pass
matplotlib.use('Qt5Agg')                

#%%
    #This one is a TOI.
ticid = 147203645
sector = 1
cam = 1
ccd = 1
path = "/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/straw-lightcurves/"
filename = "tic%012u/tic%012u_s%04u-%1u-%1u_stlc.fits" % (ticid,ticid, sector,cam,ccd)
app_bls.test_generic(filename,path)

    generate_plots_s3(ticid, sector, cam, ccd,
                  bls_bucket="/Users/smullally/TESS/lambdaSearch/test/tesssearchresults/", \
                   detrend_bucket="/Users/smullally/TESS/lambdaSearch/test/tesssearchresults/", \
                   ffilc_bucket="/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/straw-lightcurves/", \
                   outpath="/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/testing/")
except NoSuchKey:
    pass    
#%%
#EB example
ticid = 32155340
sector = 1
cam = 1
ccd = 3
path = "/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/straw-lightcurves/"
filename = "tic%012u/tic%012u_s%04u-%1u-%1u_stlc.fits" % (ticid,ticid, sector,cam,ccd)
app_bls.test_generic(filename,path)
generate_plots_s3(ticid, sector, cam, ccd,
                  bls_bucket="/Users/smullally/TESS/lambdaSearch/test/tesssearchresults/", \
                   detrend_bucket="/Users/smullally/TESS/lambdaSearch/test/tesssearchresults/", \
                   ffilc_bucket="/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/straw-lightcurves/", \
                   outpath="/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/testing/")
#%%
#sector 16
#Planet with outliers
ticid = 229760249
ticid = 233738446
sector = 16
cam = 3
ccd = 1
path = "/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/straw-lightcurves/"
filename = "tic%012u/tic%012u_s%04u-%1u-%1u_stlc.fits" % (ticid,ticid, sector,cam,ccd)
app_bls.test_generic(filename,path)
data=generate_plots_s3(ticid, sector, cam, ccd,
                  bls_bucket="/Users/smullally/TESS/lambdaSearch/test/tesssearchresults/", \
                   detrend_bucket="/Users/smullally/TESS/lambdaSearch/test/tesssearchresults/", \
                   ffilc_bucket="/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/straw-lightcurves/", \
                   outpath="/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/testing/")
#%%
ticid = 29753525
sector = 1
cam = 1
ccd = 4
path = "/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/straw-lightcurves/"
filename = "tic%012u/tic%012u_s%04u-%1u-%1u_stlc.fits" % (ticid,ticid, sector,cam,ccd)
app_bls.test_generic(filename,path)
generate_plots_s3(ticid, sector, cam, ccd,
                  bls_bucket="/Users/smullally/TESS/lambdaSearch/test/tesssearchresults/", \
                   detrend_bucket="/Users/smullally/TESS/lambdaSearch/test/tesssearchresults/", \
                   ffilc_bucket="/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/straw-lightcurves/", \
                   outpath="/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/testing/")
#%%
ticid = 2028097053
sector = 1
cam = 1
ccd = 2
path = "/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/straw-lightcurves/"
filename = "tic%012u/tic%012u_s%04u-%1u-%1u_stlc.fits" % (ticid,ticid, sector,cam,ccd)
app_bls.test_generic(filename, path)
generate_plots_s3(ticid, sector, cam, ccd,
                  bls_bucket="/Users/smullally/TESS/lambdaSearch/test/tesssearchresults/", \
                   detrend_bucket="/Users/smullally/TESS/lambdaSearch/test/tesssearchresults/", \
                   ffilc_bucket="/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/straw-lightcurves/", \
                   outpath="/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/testing/")
#%%
#lc from PL bert FFI run
ticid = 1129033
sector = 4
cam = 1
ccd = 1
path = "/Users/smullally/TESS/lambdaSearch/bertLC/"
filename = "tic%012u/tic%012u_s%04u-%1u-%1u_lcc.fits" % (ticid,ticid, sector,cam,ccd)
app_bls.test_generic(filename, path)
generate_plots_s3(ticid, sector, cam, ccd,
                  bls_bucket="/Users/smullally/TESS/lambdaSearch/test/tesssearchresults/", \
                   detrend_bucket="/Users/smullally/TESS/lambdaSearch/test/tesssearchresults/", \
                   ffilc_bucket="/Users/smullally/TESS/lambdaSearch/bertLC/", \
                   outpath="/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/testing/")
#%%
#lc from PL bert FFI run
ticid = 402026209
sector = 2
cam = 2
ccd = 4
path = "/Users/smullally/TESS/lambdaSearch/bertLC/"
filename = "tic%012u/tic%012u_s%04u-%1u-%1u_lcc.fits" % (ticid,ticid, sector,cam,ccd)
app_bls.test_generic(filename, path)
generate_plots_s3(ticid, sector, cam, ccd,
                  bls_bucket="/Users/smullally/TESS/lambdaSearch/test/tesssearchresults/", \
                   detrend_bucket="/Users/smullally/TESS/lambdaSearch/test/tesssearchresults/", \
                   ffilc_bucket="/Users/smullally/TESS/lambdaSearch/bertLC/", \
                   outpath="/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/testing/")
#%%
#Example with lots of outliers.
ticid = 2028237787
sector = 1
cam = 1
ccd = 3
generate_plots_s3(ticid, sector, cam, ccd,
                  bls_bucket="/Users/smullally/TESS/lambdaSearch/test/tesssearchresults/", \
                   detrend_bucket="/Users/smullally/TESS/lambdaSearch/test/tesssearchresults/", \
                   ffilc_bucket="/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/straw-lightcurves/", \
                   outpath="/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/testing/")

#%%
#Examples with a second noisy section
ticid = 302296676
ticid = 302296969
ticid = 302259268
ticid = 302333689
sector = 1
cam = 1
ccd = 4
generate_plots_s3(ticid, sector, cam, ccd,
                  bls_bucket="/Users/smullally/TESS/lambdaSearch/test/tesssearchresults/", \
                   detrend_bucket="/Users/smullally/TESS/lambdaSearch/test/tesssearchresults/", \
                   ffilc_bucket="/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/straw-lightcurves/", \
                   outpath="/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/testing/")
#%%
#Run all of them 
def run_all():
    runfile="/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/s0001-kdwarf/runfile.txt"
    info = np.loadtxt(runfile)
    ticids = info[:,1]
    sectors = info[:,2]
    cams= info[:,3]
    ccds = info[:,4]
    
    for j, tic in enumerate(ticids[45:]):
        try:
            generate_plots_s3(int(tic), sectors[j], cams[j], ccds[j])
            plt.close('all')
        except:
            print("Missing %i" % int(tic))
            pass


#%%
blsfile = "/Users/smullally/TESS/lambdaSearch/strawTests/blsResults/s0001-kdwarf/all_bls_results.csv"

results = p.read_csv(blsfile,comment='#')
results['id'] = list(map(lambda x: "%015.02f" % x, results['tic.pn']))
print(results.columns)

#%%
want = (results['bls_snr'] > 3) & \
        (results['n_trans_good'] > 2)
results.sort_values('n_trans_good', ascending=False, inplace=True)

print(results[want][['id','bls_period_day','bls_depth_frac','bls_snr','n_trans_good']])

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

ticid = 92560718
camera = 1
ccd = 4
sector = 1
one_result = results.loc(575)

lc_path = "tic%012u/" % ticid
search_path = "tic%012u/" % ticid
rootname = "tic%012u_s%04u-%1u-%1u" %  (ticid, sector, camera, ccd)
bls_name = "%s_plsearch.csv" % rootname
det_name = "%s_detrend.fits" % rootname
lc_name = "%s_stlc.fits" % rootname

print(bls_name,det_name,lc_name)


#%%

hdu_lc = loadFitsFromUri(lc_bucket, lc_path, lc_name)

head = hdu_lc[1].header
plt.figure()
plt.imshow(hdu_lc[2].data, vmax=400)
plt.plot(head['CUBECOL'], head['CUBEROW'],'o',mfc="None", mec='red',mew=1)
plt.figure(figsize=(14,3))
plt.plot(hdu_lc[1].data['TIME'],hdu_lc[1].data['SAP_FLUX'],'.')

#%%
hdu = loadFitsFromUri(search_bucket, search_path, det_name)
#%%
time = hdu[1].data['TIME']
flux = hdu[1].data['DETREND_FLUX']

plt.figure(figsize=(14,3))
plt.plot(time, flux,'.')

plot_bls_folds(time,flux,one_result,ticid, zoom=True)

        