#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 10:04:11 2021

@author: smullally

Data analysis script for tess false alarm project.
Goal here is to match the  known planets with those that are found for each pipeline HLSP.
"""

import pandas as p
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

results_dir = "/Users/smullally/Science/tess_false_alarms/vet_results/joe_mar_12/"
result_file = "/Users/smullally/Science/tess_false_alarms/vet_results/joe_mar_12/summary_tces.csv"

#spoc
planet_target_file = "spocFilenames_deepshort_ebplanets_mag14_v2.txt"
noplanet_target_ffile = "spocFilenames_noebplanets_mag13.txt"

#qlp
planet_target_file = "qlpFilenames_deepshort_ebplanets_mag14_v2.txt"
noplanet_target_ffile = "qlpFilenames_noebplanets_mag13.txt"

#tic for stars we care about.
tic_star_file = "/Users/smullally/Science/tess_false_alarms/keplerTargets/target_selection/tic_kic_stars_15.csv"
tic = p.read_csv(tic_star_file, sep = ',', usecols=['ID','KIC','Tmag','Vmag','Jmag','GAIA','logg','Teff'])


#%%
def get_planet_list(ddir = "/Users/smullally/Science/tess_false_alarms/keplerTargets/catalogs/",
                       max_period=25, min_depth=20):

    #%Return data frame of EBs and planets with periods less than max_period
    #planets
    koifile = ddir + "cumulative.csv"
    data = ascii.read(koifile)
    kois = data.to_pandas()
    kois["koi_time0tjd"] = kois["koi_time0bk"] + 2454833.0 - 2457000.0
    deepshortplanets = kois[(kois.koi_depth > min_depth) & (kois.koi_period <= max_period)]
    columns = ['kepid','kepoi_name', 'kepler_name', 'koi_pdisposition',
               'koi_period','koi_period_err1', 'koi_time0bk','koi_time0tjd',
               'koi_depth','koi_steff']
    
    planets = deepshortplanets[columns]
    #EBs -- ebs don't have depth, so let's ignore for themoment.
    #ebfile = ddir + "villanova_eb_catalog.csv"
    #ebs = p.read_csv(ebfile,sep=',',comment="#")
    #shortebs = ebs[(ebs.period<25)]
    
    
    return planets


#%%
#Goal is to look at targets with known planets and see if we find them again.

planets = get_planet_list()
planetstic = planets.merge(tic, how='inner',left_on='kepid',right_on='KIC')

#%%

def get_corazon_results(summary_filename):
    """
    Fill a data frame with the period, epoch depth found by corazon.
    expect a csv file made from catting together the .csv files produced
    by corazon
    
    """
    columns=['ticid', 'tceno', 'sector', 'period','epoch', 
             'depth', 'duration','snr','disp','reason']
    
    tcedf = p.read_csv(summary_filename, sep=",", names=columns)
    tcedf['epoch_tjd'] = tcedf['epoch'] + 2454833.0 - 2457000.0   
    
    
    tcedf['tce_name'] = tcedf['ticid'].astype(str) + '_'+tcedf['tceno'].astype(str) + '_'+ tcedf['sector'].astype(str)

    return tcedf

#%%

qlp_summary_file = "/Users/smullally/Science/tess_false_alarms/vet_results/joe_mar_20/second-dataset/qlp/ebplanets-mag14-v2/tce_summary_qlp.csv"
qlp_tces = get_corazon_results(qlp_summary_file)



#%%
#for each planet determine if there are results in qlp_tces, also deteremine if those results
#Match the expected period and epoch. return the difference in epoch  (worried about TTVs)

def match_signals(orig_list, tces, err = 0.1, sector=14):
    """
    Return 2 lists of arrays the length of orig_list
    first is was it searched in that sector (1/0)
    second whether it was found in that sector.
    
    err is a percentage error on the period


    """
    searched = np.zeros(len(orig_list))
    found = np.zeros(len(orig_list))
    found_name =  np.array(found, dtype=str)
    
    for i, planet in orig_list.iterrows():
        
          ticid = planet['ID']
          period = planet['koi_period']
          
          tce_stars = tces[ tces['ticid'] == 'TIC ' + str(ticid) ]
          #print(tce_stars)
          
          if sector in list(np.unique(tce_stars['sector'])):
              searched[i] = 1
          
          for j,tce in tce_stars.iterrows():
             tce_period = tce['period']
             
             if ((1-tce_period/period) < err) | ((1-2*tce_period/period)< err):
                 found[i] = 1
                 found_name[i] = tce['tce_name']
                 
    return searched, found, found_name
          
#%%
searched14,found14,found_name14 = match_signals(planetstic, qlp_tces, err=0.03,sector=14)
searched15,found15,found_name15 = match_signals(planetstic, qlp_tces, err=0.03,sector=15)
searched16,found16,found_name16 = match_signals(planetstic, qlp_tces, err=0.03,sector=16)
searched26,found26,found_name26 = match_signals(planetstic, qlp_tces, err=0.03,sector=26)

#%%
import lightkurve as lk

allsearched = searched14 + searched15 + searched16 + searched26
keptic_searched = planetstic[allsearched>0]
keptic_notsearched = planetstic[allsearched==0]
allfound = found14 + found15+ found16 + found26
keptic_searched_found = planetstic[(allsearched>0) & (allfound>0)]
keptic_searched_notfound = planetstic[(allsearched>0) & (allfound==0)]

for i, target in keptic_searched_found[(keptic_searched_found['koi_depth']<2000) & \
                                          (keptic_searched_found['koi_period']<10)].iterrows():
    
    name = 'TIC ' + str(target['ID'])
    print(target[['ID','kepoi_name','koi_pdisposition','koi_period','koi_depth','Tmag']])
    tcefound = qlp_tces[qlp_tces['ticid'] == name]
    print(tcefound[['ticid','sector','period','depth','snr','disp','reason']])
    sectors = np.unique(tcefound['sector'])
    plt.close()
    for s in sectors:

        lk.search_lightcurve(name, mission='TESS', sector=s, author="qlp")[0].download().remove_nans().plot()
        plt.pause(.05)
        
    input("Press Enter to continue...")

          
#%%
columns = ['kepid','kepoi_name', 'koi_period','koi_pdisposition', 'koi_time0tjd','koi_depth',
               'ticid', 'tceno', 'sector', 'period','epoch', 
                'depth', 'duration','snr','disp','reason',
                'found', 'depoch','dperiod']
tce_columns=['ticid', 'tceno', 'sector', 'period','epoch', 
             'depth', 'duration','snr','disp','reason']
    