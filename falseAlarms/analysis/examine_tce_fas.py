#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 11:48:53 2021

@author: smullally

This code looks at the qlp TCEs we made and plots up a cumulative frequency
above a certain snr.
It also pulls in TIC information including contamination fraction to investigate
options. Contratio is from the TIC where some of the nans are calculated by 
counting the number stars within 1 aperture radii and 1.5magnitude difference.

"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as p
import scipy.stats as stats

qlp_tces = "/Users/smullally/Science/tess_false_alarms/vet_results/joe_mar_2021/first-dataset/qlp_noebplanets_tcesum.csv"
spoc_tces  = "/Users/smullally/Science/tess_false_alarms/vet_results/joe_mar_2021/first-dataset/spoc_noebplanets_tcesum.csv"

col_names = ["tic","pn","sector","period","epoch","depth","dur","snr","disp","reason"]
qlp_df = p.read_csv(qlp_tces, names=col_names)
spoc_df = p.read_csv(spoc_tces, names=col_names)

qlp_df['uniqueid'] = qlp_df['tic'].str[4:].astype(int)*100+qlp_df['sector']

#Created using uncrowded_kepler_targets.py
targetdf = p.read_csv("/Users/smullally/Science/tess_false_alarms/keplerTargets/target_selection/target_tic_contamination_20210528.txt")
targetdf['uniqueid'] = targetdf['ticid']*100 + targetdf['sector']

#Merge target information into qlp TCEs.
qlp_tces = p.merge(qlp_df, targetdf[['Tmag', 'Hmag', 'Vmag','contratio','aperture','uniqueid']], left_on="uniqueid", right_on="uniqueid", how="left" )
#%%
#Plot distributions of the snr of the TCEs.

bins = np.arange(1,100,.2)

plt.figure()
n,bins, patches = plt.hist(qlp_tces['snr'],bins = bins,histtype='step')
plt.yscale('log')
#Passing
passed = qlp_tces['disp'] == 'PASS'
failing = qlp_tces['disp'] == 'FAIL'

plt.figure()
n, bins, patches = plt.hist(qlp_tces['contratio'], histtype='step', bins=np.arange(0,3,.01))
plt.title('Contrast Ratio')


uncrowded = qlp_tces['contratio'] <= 0.1
bright = qlp_tces['Tmag'] < 11

plt.figure()
n, bins, patchees = plt.hist(qlp_tces[uncrowded]['snr'], histtype='step',bins=np.arange(1,100,.2))
plt.title('uncrowded snr')
plt.yscale('log')

#%%
#Create a data frame o just the first TCE found in each
#Look at the cumulative distribution of those

def plot_cumfreq(df, over = 10):
    
    maxsnr = int(np.ceil(np.max(df['snr'])))
    plt.figure(figsize=(6,8))
    plt.subplot(211)
    bins = np.arange(0,maxsnr,1/over)
    n, abins, patches = plt.hist(df['snr'], bins = bins, histtype='step')
    plt.yscale('log')
    plt.xlim(0,20)
    plt.title('First BLS Detection Histogram')
    
    cumres = stats.cumfreq(df['snr'], maxsnr*over, (0, maxsnr))
    
    xbins = cumres.lowerlimit + np.linspace(0, cumres.binsize*cumres.cumcount.size,
                                           cumres.cumcount.size)
    cumfract = 100*(cumres.cumcount[-1]-cumres.cumcount)/len(df)
    plt.subplot(212)
    plt.step(xbins, cumfract)
    plt.yscale('log')
    plt.xlim(0,20)
    plt.title('First BLS Detection Cumulative Frequency')
    plt.ylabel('Fraction Found with Greater SNR')
    
    nearestx = xbins[np.argmin(np.abs(cumfract-0.3))]
    print(nearestx)
    plt.hlines(0.3, 0, nearestx, color='red')
    plt.annotate(str(nearestx), (nearestx, 0.3), color='red')

    return(cumfract)

#%%
#Just the first tce found in the search
qlp_tces_first = qlp_tces.query('pn == 1').copy()

cumfract = plot_cumfreq(qlp_tces_first)
plt.title("All QLP First Cum Distribution")


cumfract = plot_cumfreq(qlp_tces_first[uncrowded])
plt.title("Uncrowded QLP First Cum Distribution")

#%
#Look at differences based on period.

cumfract = plot_cumfreq(qlp_tces_first.query('period > 8'))
plt.title('Period > 8d)')

cumfract = plot_cumfreq(qlp_tces_first.query('period < 4'))
plt.title('Period < 4d)')

cumfract = plot_cumfreq(qlp_tces_first[uncrowded].query('period > 2 & period < 9'))
plt.title('Uncrowded and Periods between 2 and 9')

#%%
#Pull in simulated TCEs.
#See what that looks like.
sim_datafiles = ["/Users/smullally/Science/tess_false_alarms/simulations/simulated_tces_200_1200_20210215_v1.txt",
                 "/Users/smullally/Science/tess_false_alarms/simulations/simulated_tces_200_1200_20210215_v2.txt",
                 "/Users/smullally/Science/tess_false_alarms/simulations/simulated_tces_200_400_20210219_v3.txt"]

names = ['ID', 'pn', 'sector', 'period', 'epoch', 'depth', 'duration', 'snr', 'disp', 'reason', 'noise']

sim = p.DataFrame(columns=names)
for file in sim_datafiles:
    simdf = p.read_csv(file, names=names, skiprows=1)
    sim = sim.append(simdf)

cumfract = plot_cumfreq(sim)
plt.title('Gaussian noise')


#%%

#Vet the high snr ones to throw out obvious EBs.
#These ones are uncrowded, so have no reason to show a TCE.
#also high snr.
#21 in total
qlp_high_snr = qlp_tces_first[uncrowded & (qlp_tces_first.snr > 5)].copy()

qlp_ois = qlp_tces_first[uncrowded].query('period > 2 & period < 9 & snr > 5')
#Calculate the number of transits found

#%%
import exovetter.vetters as vet

#Vetting Instructions -- But need ot import other thing before we can run this.
i = 7;
alc, atpf = get_lk(qlp_ois.iloc[i]); 
atce = make_tce(qlp_ois.iloc[i], offset= const.btjd)
cent.run(atce, atpf, plot=True)
fold.run(atce, alc, plot=True)
qlp_ois.iloc[i]