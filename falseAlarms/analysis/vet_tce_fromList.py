#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 17:18:09 2021

@author: smullally
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as p
import exovetter
from exovetter.tce import Tce
from exovetter import const 
import astropy.units as u
import lightkurve as lk
import corazon as crz
import exovetter.vetters as vet

def make_tce(adf, offset = 0 *u.day):
    
    atce = Tce(period = adf['period'] * u.day,
               epoch = adf['epoch'] *u.day,
               epoch_offset = offset,
               duration=adf['dur'] * u.hr, 
               depth=adf['depth'] * const.ppm,
               snr = adf['snr'])
    
    return atce


def vet_tce(tce, lc, tpf):
    """Pull up plots of the full and folded light curve.
        There are all plot from exovetter
    """
    
    tc = vet.TransitPhaseCoverage()
    tpc = tc.run(tce,lc)
    
    try:
        mod = vet.ModShift()
        modshift = mod.run(tce, lc)
        modshift.plot()
    except:
        pass
    
    sweet = vet.Sweet()
    sweetvet = sweet.run(tce,lc)
    sweet.plot()
    
    cent = vet.Centroid()
    centout = cent.run(tce, tpf, plot=True)
    
    oe = vet.OddEven()
    oddeven = oe.run(tce,lc)
    oe.plot()
    
    fold = vet.VizTransits(smooth=8, max_transits=5)
    ntransits = fold.run(tce, lc, plot=True)
    
    return oddeven, ntransits, centout

def get_lk(tcedf, author = "qlp", mission = "TESS", size = 11):
    """Returns lc and tpf given a dataframe line"""
    
    ticid = tcedf['tic'][4:]
    sector = tcedf['sector']
    
    lc = crz.gen_lightcurve.hlsp(ticid, sector, author=author)
    tpf = lk.search_tesscut(f"TIC {ticid}", sector = sector).download(cutout_size = size)
    
    return lc, tpf
#%%

#Get list of TCEs we could vet from the summary file found in the top level directory
qlp_tces = "/Users/smullally/Science/tess_false_alarms/vet_results/joe_mar_2021/first-dataset/qlp_noebplanets_tcesum.csv"
col_names = ["tic","pn","sector","period","epoch","depth","dur","snr","disp","reason"]
qlp_df = p.read_csv(qlp_tces, names=col_names)
#%%
#Get tces for df
qlp_tces = "/Users/smullally/Science/tess_false_alarms/vet_results/joe_mar_2021/second-dataset/qlp_ebplanets_tcesum.csv"
col_names = ["tic","pn","sector","period","epoch","depth","dur","snr","disp","reason"]
qlpdf = p.read_csv(qlp_tces, names=col_names)


#%%
#Function so that I only need to specify theline number of the file for maual vetting.

def vet_number(n, tcedf, update = None):
    
    num= n
    #ticid = tcedf.iloc[num].tic
    #sector = tcedf.iloc[num].sector
    atcedf = tcedf.iloc[num]
    alc, atpf = get_lk(atcedf)
    
    atce = make_tce(atcedf, offset= const.btjd)
    
    if update is not None:
        for key in update:
            atce[key] = update[key]
    
    print(atce)
    
    oddeven, ntransits, centout = vet_tce(atce, alc, atpf)
    
    alc.plot()
    
    print(atce)
    print(ntransits)
    
    return alc, atce, atpf



   

#%%
want = qlp_df['period'] < 10
tcedf = qlp_df[want].query('snr >8')
#%% 
n = 0
print(tcedf.iloc[n][['tic','pn','disp','depth','snr','reason']] ) 
alc, atce, atpf = vet_number(n,tcedf, update=None)

#%%
mod = vet.ModShift()
modshift = mod.run(atce, alc, plot=True)


#%%
#Other examples..

ticid = "TIC 120576881"
sector = 14
tcedf = qlp_df.query('tic == @ticid & sector == @sector')
print(tcedf[['tic','pn','sector','period','snr','disp','reason']])
 #%%
choice = 0
alc, atpf = get_lk(tcedf.iloc[choice])
atce = make_tce(tcedf.iloc[choice], offset = const.btjd)
vet_tce(atce, alc, atpf)