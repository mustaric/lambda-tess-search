#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 14:30:18 2022

@author: smullally
"""

#Need to determine if the ephemeris matches a known KOI or EB.

#Table
from astropy.io import ascii
from astropy.table.table import Table
import numpy as np
from astroquery.mast import Catalogs
from exovetter.tce import Tce
from exovetter import const as exo_const
from astropy import units as u



#%%
#Reformat to a standard set of information.
#So easily access period and epoch for both catalogs

    
def reformat_variable_catalogs():
    koi_file = "/Users/smullally/Science/tess_false_alarms/keplerTargets/catalogs/cumulative.csv"
    eb_file = "/Users/smullally/Science/tess_false_alarms/keplerTargets/catalogs/villanova_eb_catalog.csv"
    
    kois = ascii.read(koi_file)
    #koi_tces = kois['kepid','kepoi_name', 'koi_period', 'koi_time0bk', 'koi_depth']
    
    ebs = ascii.read(eb_file, header_start=0, delimiter=',')
    
    star = list(kois['kepid']) + list(ebs['KIC'])
    signal = list(kois['kepoi_name']) + list(ebs['KIC'].astype(str))
    period = list(kois['koi_period']) + list(ebs['period'])
    epoch = list(kois['koi_time0bk']) + list(ebs['bjd0']-54833.0)
    catalog = list(['koi'] * len(kois)) + list(['eb'] * len(ebs))
    
    vartable = Table([star, signal, period, epoch, catalog],
                    dtype = [int, str, float, float, str],
                    names = ['KIC','signal','period_day','epoch_bkjd','catalog'])
    varname = "/Users/smullally/Science/tess_false_alarms/keplerTargets/catalogs/varTable.csv"
    vartable.write(varname, format='csv', overwrite=False)

#%%
#For a list of KIC Ids return known EBs or planet IDs from a catalog
#vartab is reading in the table created above (vartable)
def get_tces(kepid_list, vartab):
    
    tces = list()
    #type(tces)
    #print(kepid_list)
    for KIC in kepid_list:
        index = np.where(vartab['KIC'] == int(KIC))
        for idx in index[0]:
            atce = Tce(period = vartab[idx]['period_day'] * u.day, 
                       epoch = vartab[idx]['epoch_bkjd'] * u.day,
                       epoch_offset = exo_const.bkjd,
                       depth = 1000 * exo_const.ppm,
                       duration = 5 * u.hr,
                       kicid = KIC,
                       signal_name = vartab[idx]['signal'])

            tces.append(atce)

    return tces
     

#%%

#For a TIC id.  Find all nearby KIC Ids  (within ~1arc min)

def nearby_kic(name, sep_deg=90/3600):
    
    #name = "TIC %u" % ticid
    catalog_data = Catalogs.query_object(name, radius = sep_deg, catalog="TIC")
    #print(catalog_data.columns)
    #print(catalog_data['ID','Tmag','Teff','logg','KIC'])
    
    return catalog_data['ID','Tmag','Teff','logg','KIC']       

#%%
# Determine if any of the Kepler tce events occur at the same time as the TESS one.
#

def transit_time_match(periods, epochs, time_range, time_error):
    #For two period/epoch pairs
    #Determine if the transit times overlap within the time_error given the time_range.
    N1 = np.floor((time_range[0]-epochs[0]) / periods[0])
    N2 = np.ceil((time_range[1]-epochs[0]) / periods[0])
    transit_times1 = np.arange(N1, N2, 1) * periods[0] + epochs[0]
    
    
    N1 = np.floor((time_range[0]-epochs[1]) / periods[1])
    N2 = np.ceil((time_range[1]-epochs[1]) / periods[1])
    transit_times2 = np.arange(N1, N2, 1) * periods[1] + epochs[1]
    
    all_diffs = list()
    
    for tt1 in transit_times1:
        for tt2 in transit_times2:
            all_diffs.append(np.abs(tt1.value-tt2.value))

    all_diff_array = np.array(all_diffs)
    small_diffs = all_diff_array[all_diff_array <= time_error.to(time_range[0].unit).value]
    num_too_close = len(small_diffs)
    
    return(num_too_close, small_diffs)

#%%
#Main code runs the above.
#Goal is for each to ID as, 
#1) star + transit time match known = PC, 
#2) transit time matches nearby star = CO
#3) does not match a known ephemeris
# For OC ad PC we need to include the star and period and depth.
#

def get_true_nature(atce):
    """
    unkown means does not match anything in the varible files provided.
    """
    
    known_variables_file = "/Users/smullally/Science/tess_false_alarms/keplerTargets/catalogs/varTable.csv"
    vartab = ascii.read(known_variables_file)
    
    nearby_stars = nearby_kic(atce['star'])
    nearstars = nearby_stars[~nearby_stars['KIC'].mask]

    closekic = nearstars[0]['KIC']
    
    #Get all TCEs associated with the nearby stars    
    nearby_tces = get_tces(nearstars['KIC'], vartab)
    
    status = ""
    
    for neartce in nearby_tces:
        
        periods = [atce['period'], neartce['period']]
        epochs = [atce['epoch'] - atce['epoch_offset'], 
                  neartce['epoch'] - neartce['epoch_offset']]
        
        time_range = [1711*u.d - exo_const.btjd, 1763*u.d - exo_const.btjd]
        time_error = atce['duration']*0.5
        
        num, small_diffs = transit_time_match(periods, epochs, time_range, time_error)

        if num>=1:
            if closekic == neartce['kicid']:
                status = ';'.join([status, "PC: %s, num:%u" % (neartce['signal_name'], num)])
            else:
                status = ';'.join([status, "CO: %s, num:%u" % (neartce['signal_name'], num)])

    #import pdb;pdb.set_trace()
    if status == "":
        status = "unknown"
    
    return(status)
    

#%% 
    
from exovetter import utils
def example():
    atoi = utils.get_mast_tce('K00752.01')
    koi752 = atoi[0]
    koi752['star'] = 'KIC 10797460'
    print(koi752)
    status = get_true_nature(koi752)
    
    print(status)
    
    atoi = utils.get_mast_tce('K00759.01')
    koi759 = atoi[0]
    koi759['star'] = 'KIC 11018648'
    print(koi759)
    status = get_true_nature(koi759)
    print(status)
    
    atoi = utils.get_mast_tce('K00759.01')
    koi759 = atoi[0]
    koi759['star'] = 'KIC 11018636'  #fake for a CO test
    print(koi759)
    status = get_true_nature(koi759)
    print(status)
    
    