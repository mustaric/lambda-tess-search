#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 16:38:18 2022

@author: smullally
"""
import sys
sys.path.append("/Users/smullally/Python_Code/parmap/")
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import parmap  #This code lets me speed things up by doing archive callsin parallel
import warnings
warnings.filterwarnings('ignore')

#%%

#Useful Functions
import exovetter.vetters as vet
import corazon as crz
import lightkurve as lk
from exovetter.tce import Tce
from exovetter import const 
import astropy.units as u
import corazon

def get_lk(tcedf, author = "qlp", mission = "TESS", size = 11):
    """Returns lc and tpf given a dataframe line"""

    ticid = tcedf['tic'][4:]
    sector = tcedf['sector']
    
    lc = crz.gen_lightcurve.hlsp(ticid, sector, author=author)
    
    clean_lc = crz.planetSearch.clean_timeseries(lc['time'].value, 
                    lc['flux'].value, lc['quality'].value, 95, 19, 4.5, sector) 
    clean = lk.LightCurve(clean_lc[0],clean_lc[1])
    clean.meta = lc.meta
    clean.time.format = 'btjd'
    
    
    tpf = lk.search_tesscut(f"TIC {ticid}", sector = sector).download(cutout_size = size)
    
    return clean, tpf

def make_tce(adf, offset = 0 * u.day):
    
    atce = Tce(period = adf['period'] * u.day,
               epoch = adf['epoch'] *u.day,
               epoch_offset = offset,
               duration=adf['dur'] * u.hr, 
               depth=adf['depth'] * const.ppm,
               snr = adf['snr'])
    
    return atce
#%%
def vet_tce(tce, lc, tpf, plot=False):
    """Pull up plots of the full and folded light curve.
        There are all plot from exovetter
    """
    
    results = []
    
    tc = vet.TransitPhaseCoverage()
    tpc = tc.run(tce,lc)
    results.append(tpc)
    
    lpp = vet.Lpp()
    lppvet = lpp.run(tce,lc)
    results.append(lppvet)
    
    try:
        mod = vet.ModShift()
        modshift = mod.run(tce, lc)
        if plot:
            modshift.plot()
        results.append(modshift)
    except:
        pass
    
    sweet = vet.Sweet()
    sweetvet = sweet.run(tce,lc)
    if plot:
        sweet.plot()
    results.append(sweetvet)
    cent = vet.Centroid()
    centout = cent.run(tce, tpf, plot=plot)
    results.append(centout)
    
    #oe = vet.OddEven()
    #oddeven = oe.run(tce,lc)
    #oe.plot()
    #results.append(oddeven)
    
    fold = vet.VizTransits(smooth=8, max_transits=5, transit_only=True)
    ntransits = fold.run(tce, lc, plot=plot)
    results.append(ntransits)
    
    return results
    
def vet_get_results(tcedfrow, plot=False, offset=const.btjd):
    tcedf = tcedfrow[1]
    atce = make_tce(tcedf, offset=offset)
    
    try:
        alc, atpf = get_lk(tcedf) 
    except:
        return dict(tcedf)
    
    if atpf is not None:
        try:
            results = vet_tce(atce, alc, atpf, plot=plot)
        except AssertionError:
            results = None
    else:
        results = None
        
    tcevet = dict(tcedf)
    if results is not None:
        for r in results:
            for k in r.keys():
                if k[0:3] != "plo":  #Don't include plot data
                    tcevet[k] = r[k]
        
    return tcevet    
#%%

def show_tce(tcedfrow, plot=False, offset=const.btjd):
    tcedf = tcedfrow
    alc, atpf = get_lk(tcedf) 
    atce = make_tce(tcedf, offset=offset)

    fold = vet.VizTransits(smooth=8, max_transits=5, transit_only=True)
    results = fold.run(atce, alc, plot=True)        
    fold = vet.VizTransits(smooth=5, max_transits=5, transit_only=False)
    results = fold.run(atce, alc, plot=plot)
    
    print(tcedf)
        
    return results  

def get_disp(match_str):
    
    if "PC" in match_str: 
        if "num:1" not in match_str:
            return "PC"
    elif "CO" in match_str:
        return "CO"
    else:
        return "FP"


#%%


notransit_file = "/Users/smullally/Science/tess_false_alarms/vet_results/joe_mar_2021/first-dataset/qlp_noebplanets_tcesum.csv.em.csv"
planet_file = "/Users/smullally/Science/tess_false_alarms/vet_results/joe_mar_2021/second-dataset/qlp_ebplanets_tcesum.csv.em.csv"
col_names = ["index","tic","pn","sector","period","epoch","depth","dur","snr","disp","reason", "ephmatch"]
notransit_df = pd.read_csv(notransit_file, names=col_names, skiprows=1)
planet_df = pd.read_csv(planet_file, names=col_names, skiprows=1)

#Created using uncrowded_kepler_targets.py
targetdf = pd.read_csv("/Users/smullally/Science/tess_false_alarms/keplerTargets/target_selection/target_tic_contamination_20210528.txt")
print(len(targetdf))
targetdf = targetdf.drop_duplicates(subset=['ticid','sector', 'Tmag'])
print(len(targetdf))

#Uniqueid
notransit_df['uniqueid'] = list(map( lambda x, y : "%s-s%u" % (x, y), notransit_df['tic'], notransit_df['sector']))
planet_df['uniqueid'] = list(map( lambda x, y : "%s-s%u" % (x, y), planet_df['tic'], planet_df['sector']))
targetdf['uniqueid'] = list(map( lambda x, y : "TIC %s-s%u" % (x, y), targetdf['ticid'], targetdf['sector']))

#Merge target information into qlp TCEs.
notransit_tces = pd.merge(notransit_df, targetdf[['Tmag', 'Hmag', 'Vmag','contratio','aperture','uniqueid']], left_on="uniqueid", right_on="uniqueid", how="left" )
planet_tces = pd.merge(planet_df, targetdf[['Tmag', 'Hmag', 'Vmag','contratio','aperture','uniqueid']], left_on="uniqueid", right_on="uniqueid", how="left" )

planet_tces['matchdisp'] = list(map(lambda x : get_disp(x), planet_tces['ephmatch']))
notransit_tces['matchdisp'] = list(map(lambda x : get_disp(x), notransit_tces['ephmatch']))


pcs = planet_tces['matchdisp'] == 'PC'
cos = notransit_tces['matchdisp'] == 'CO'
fps = notransit_tces['matchdisp'] == 'FP'
plbright = planet_tces['Tmag'] < 12
ntbright = notransit_tces['Tmag'] < 12
pluncrowd = planet_tces['contratio'] < 0.3
ntuncrowd = notransit_tces['contratio'] < 0.3

#%%
#This selection of the noise has those that mach a CO or PC removed. (unless only matches 1 transit)

bestnoise = notransit_tces[fps & ntbright]

numunique = len(np.unique(bestnoise['uniqueid']))
print(numunique)


plt.figure()
bins = np.arange(1,70,.5)
n,bins, patches = plt.hist(bestnoise['snr'], bins= bins, histtype='step')
plt.yscale('log')
plt.xlabel('SNR')
plt.title('bright')


bestnoise = notransit_tces[fps & ntbright & ntuncrowd]

plt.figure()
bins = np.arange(1,70,.5)
n,bins, patches = plt.hist(bestnoise['snr'], bins= bins, histtype='step')
plt.yscale('log')
plt.xlabel('SNR')
plt.title('bright and uncrowded')

#%%
#Let's do a quick vet on all of these so that we can apply an lpp cut.

for b in np.arange(2, 5,.1):
    snrrange = (notransit_tces['snr']>b) & (notransit_tces['snr'] <= (b+.1))
    bestnoise = notransit_tces[snrrange & fps & ntbright]
    print(len(bestnoise))
    print("b is %f" % b)
    result = parmap.parmap(vet_get_results, bestnoise.iterrows(),engine='serial')
    new = list(filter(lambda x : x is not None, result))
    vetted_notransit = pd.DataFrame(new)
    vetted_notransit.to_pickle("/Users/smullally/Science/tess_false_alarms/notransit_vets/notransit_highsnr_fps_20220513_%u.pickle" % (int(b*10)))


