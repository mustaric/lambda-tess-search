#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The goal of this code is to compare the vetting thersholds for 
both the planet and the 

Created on Wed May 11 09:52:50 2022

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

def get_lk(tcedf, author = "qlp", mission = "TESS", size = 11):
    """Returns lc and tpf given a dataframe line"""

    ticid = tcedf['tic'][4:]
    sector = tcedf['sector']
    
    lc = crz.gen_lightcurve.hlsp(ticid, sector, author=author)
    
    clean_lc = crz.planetSearch.clean_timeseries(lc['time'].value, lc['flux'].value, lc['quality'].value, 95, 19, 4.5, sector) 
    clean = lk.LightCurve(clean_lc[0],clean_lc[1])
    clean.meta = lc.meta
    clean.time.format = 'btjd'
    
    #tpf = lk.search_tesscut(f"TIC {ticid}", sector = sector).download(cutout_size = size)
    
    return clean
"""
def get_lk(tcedf, author = "qlp", mission = "TESS", size = 11):
    #Returns lc and tpf given a dataframe line

    ticid = tcedf['tic'][4:]
    sector = tcedf['sector']
    
    lc = crz.gen_lightcurve.hlsp(ticid, sector, author=author)
    
    tpf = lk.search_tesscut(f"TIC {ticid}", sector = sector).download(cutout_size = size)
    
    return lc, tpf
"""

def make_tce(adf, offset = 0 * u.day):
    
    atce = Tce(period = adf['period'] * u.day,
               epoch = adf['epoch'] *u.day,
               epoch_offset = offset,
               duration=adf['dur'] * u.hr, 
               depth=adf['depth'] * const.ppm,
               snr = adf['snr'])
    
    return atce

def vet_tce(tce, lc, plot=False):
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
    #cent = vet.Centroid()
    #centout = cent.run(tce, tpf, plot=plot)
    #results.append(centout)
    
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
    alc = get_lk(tcedf) 
    atce = make_tce(tcedf, offset=offset)
    results = vet_tce(atce, alc, plot=plot)
        
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

#%%

#Using qlp
notransit_file = "/Users/smullally/Science/tess_false_alarms/vet_results/joe_mar_2021/first-dataset/qlp_noebplanets_tcesum.csv"
planet_file = "/Users/smullally/Science/tess_false_alarms/vet_results/joe_mar_2021/second-dataset/qlp_ebplanets_tcesum.csv"
col_names = ["tic","pn","sector","period","epoch","depth","dur","snr","disp","reason", "match"]
notransit_df = pd.read_csv(notransit_file, names=col_names)
planet_df = pd.read_csv(planet_file, names=col_names)

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
notransit_tces = pd.merge(notransit_df, targetdf[['Tmag', 'Hmag', 'Vmag','contratio','aperture','uniqueid']], left_on="uniqueid", right_on="uniqueid", how="inner" )
planet_tces = pd.merge(planet_df, targetdf[['Tmag', 'Hmag', 'Vmag','contratio','aperture','uniqueid']], left_on="uniqueid", right_on="uniqueid", how="inner" )

print(len(notransit_df), len(notransit_tces), len(targetdf))
print(len(np.unique(targetdf['uniqueid'])))
print(len(planet_df), len(planet_tces))
print(len(np.unique(notransit_df['uniqueid'])), len(np.unique(notransit_tces['uniqueid'])))

#%%
#Get N tces at random for testing
#This cell can take some time depending on the size.
r = np.random.randint(0, len(notransit_tces), size = 200)
notransit_samp = notransit_tces.iloc[r]
result = parmap.parmap(vet_get_results, notransit_samp.iterrows(),engine='multi')
new = list(filter(lambda x : x is not None, result))
vetted_notransit = pd.DataFrame(new)

#%
vetted_notransit.to_pickle("/Users/smullally/Science/tess_false_alarms/notransit_df_20220520.pickle")

#%%
r = np.random.randint(0, len(planet_tces), size = 200)
planet_samp = planet_tces.iloc[r]
result = parmap.parmap(vet_get_results, planet_samp.iterrows(),engine='multi', timeout_sec=255)
bnew = list(filter(lambda x : x is not None, result))
vetted_planets = pd.DataFrame(bnew)
vetted_planets.to_pickle("/Users/smullally/Science/tess_false_alarms/planets_df_20220520.pickle")

#%%
#Read in those pickle files and investiate poopulations

vetted_planets = pd.read_pickle("/Users/smullally/Science/tess_false_alarms/planets_df_20220520.pickle")
vetted_notransit = pd.read_pickle("/Users/smullally/Science/tess_false_alarms/notransit_df_20220520.pickle")

#%%
plt.plot(vetted_planets['period'], vetted_planets['norm_lpp'], 'r.')
plt.plot(vetted_notransit['period'], vetted_notransit['norm_lpp'], 'k.')


plt.figure()
bins= np.arange(.5,20,.5)
plt.hist(vetted_planets['norm_lpp'], histtype='step', bins=bins)
plt.hist(vetted_notransit['norm_lpp'], histtype='step', bins=bins)


plt.figure()
first= vetted_planets['pn'] == 1
plt.plot(vetted_planets['pn'], vetted_planets['norm_lpp'], 'o')
#%%
plt.figure()
plt.plot(vetted_notransit['significance'], vetted_notransit['norm_lpp'],'x')
plt.plot(vetted_planets['significance'], vetted_planets['norm_lpp'],'ro')



#%%

def get_disp(match_str):
    
    if "PC" in match_str: 
        if "num:1" not in match_str:
            return "PC"
    elif "CO" in match_str:
        return "CO"
    else:
        return "FP"


#%%
# Try again, but this time use the ephemeris match inforamtion to decise which are real

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

#%%
#Define some qualifiers for the planet and notransit TCE dataframes
pcs = planet_tces['matchdisp'] == 'PC'
cos = notransit_tces['matchdisp'] == 'CO'
fps = notransit_tces['matchdisp'] == 'FP'
plbright = planet_tces['Tmag'] < 12
ntbright = notransit_tces['Tmag'] < 12
pluncrowd = planet_tces['contratio'] < 0.3
ntuncrowd = notransit_tces['contratio'] < 0.3


highsnr = planet_tces['snr'] > 2
pnone = planet_tces['pn'] == 1
bestplanets = planet_tces[pcs & pnone]
bestnoise = notransit_tces[fps & ntbright]

plt.figure()
bins = np.arange(0,15,.1)
plt.hist(bestnoise['snr'], bins=bins,histtype='step', label="noise")
plt.title("No ephemeris match snr histogram around bright stars (Tmag<12)")
print(len(bestnoise))
plt.ylabel("snr")
plt.legend()
plt.figure()
plt.hist(bestplanets['snr'], bins=bins, histtype='step', label='planets')
plt.legend()
plt.title("Planet SNR")
plt.ylabel("snr")

#%%
#Let's vet a sample of each.
bestnoise = notransit_tces[fps & ntbright]
r = np.random.randint(0, len(bestnoise), size = 400)
notransit_samp = bestnoise.iloc[r]
result = parmap.parmap(vet_get_results, notransit_samp.iterrows(),engine='multi')
new = list(filter(lambda x : x is not None, result))
vetted_notransit = pd.DataFrame(new)

#%
vetted_notransit.to_pickle("/Users/smullally/Science/tess_false_alarms/notransit_noise_20220513.pickle")

#%%
#Vet a bunch of the good planets
highsnr = planet_tces['snr'] > 1.2
pnone = planet_tces['pn'] == 1
bestplanets = planet_tces[pcs & highsnr & pnone]
print(len(bestplanets))
r = np.random.randint(0, len(bestplanets), size = 400)
planet_samp = bestplanets.iloc[r]
result = parmap.parmap(vet_get_results, planet_samp.iterrows(),engine='multi', timeout_sec=255)
bnew = list(filter(lambda x : x is not None, result))
vetted_planets = pd.DataFrame(bnew)
vetted_planets.to_pickle("/Users/smullally/Science/tess_false_alarms/planets_pc_20220513.pickle")

#%%
#Vet a bunch of centroid offset ones
pnone = notransit_tces['pn'] == 1
bestco = notransit_tces[cos & pnone]
r = np.random.randint(0, len(bestco), size = 300)
co_samp = bestco.iloc[r]
result = parmap.parmap(vet_get_results, co_samp.iterrows(), engine='multi')
new = list(filter(lambda x : x is not None, result))
vetted_co = pd.DataFrame(new)
vetted_co.to_pickle("/Users/smullally/Science/tess_false_alarms/notransit_co_20220513.pickle")

#%%
#Read in those pickle files and investiate poopulations

vetted_planets = pd.read_pickle("/Users/smullally/Science/tess_false_alarms/planets_pc_20220513.pickle")
vetted_notransit = pd.read_pickle("/Users/smullally/Science/tess_false_alarms/notransit_noise_20220513.pickle")
vetted_co = pd.read_pickle("/Users/smullally/Science/tess_false_alarms/notransit_co_20220513.pickle")
plt.figure()
bins = np.arange(0,.1, .001)
plt.hist(vetted_planets['raw_lpp'],histtype='step',label='planets', bins=bins)
plt.hist(vetted_notransit['raw_lpp'], histtype='step', label='notransit', bins=bins)
plt.legend()
plt.title("raw lpp")

plt.figure()
bins = np.arange(0,25, .2)
plt.hist(vetted_planets['norm_lpp'],histtype='step',label='planets', bins=bins)
plt.hist(vetted_notransit['norm_lpp'], histtype='step', label='notransit', bins=bins)
plt.legend()
plt.title('norm lpp')

#Conclusion:  Planets fall below norm_lpp=4.5  A few EBs will be lost this way.
#Note all those I inspect with norm_lpp > 5 are EBs, usually contact EBs. 
#Frequently found at the wrong period.

#%%
#centroid offset
cosig = list( map ( lambda x : float(x), vetted_co['significance']))
plsig = list( map ( lambda x : float(x), vetted_planets['significance']))

plt.figure()
bins=np.arange(0,1,.02)
plt.hist(plsig,histtype='step',label='planets', bins=bins)
plt.hist(cosig,histtype='step',label='co', bins=bins)
plt.legend()

#Conclsion is that sometimes there is an extra star that pulls the out of transit
#Can't just do a simple cut on significance.

#%%
#modshift -- complicated.
#This one is a bit like lpp to get rid of junk.

plt.figure()

plt.hist(vetted_planets[''])









