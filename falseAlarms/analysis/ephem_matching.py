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

def nearby_kic(name, sep_deg=63/3600):
    
    #name = "TIC %u" % ticid
    catalog_data = Catalogs.query_object(name, radius = sep_deg, catalog="TIC")
    #print(catalog_data.columns)
    #print(catalog_data['ID','Tmag','Teff','logg','KIC'])
    
    return catalog_data['ID','Tmag','Teff','logg','KIC']

def nearby_kic_local(ticid, neardir):
    
    filename = neardir + str(ticid)[:3]+"/"+"tic"+str(ticid)+"_nearstars.csv"
    df = p.read_csv(filename)
    
    data = Table.from_pandas(df)
    
    return data
    
    

#%%
#Create a nearby KIC file for each target.  This can be run using parmap.py
import os
def nearby_kic_make_local(outdir, name):
    sep_deg=90/3600
    cat = nearby_kic("TIC" + str(name), sep_deg=sep_deg)
    print(outdir + str(name)[:3]+"/"+"tic"+str(name)+"_nearstars.csv")
    if not os.path.exists(outdir + str(name)[:3]):
        os.mkdir(outdir + str(name)[:3])
        
    cat.to_pandas().to_csv(outdir + str(name)[:3]+ "/" + "tic"+ str(name) + "_nearstars.csv")
    
    if len(cat)>1:
        return 1 
    else:
        return 0

#Run this once to get out all the information needed.
import parmap
def make_local_nearby():
    outdir = list(["/Users/smullally/Science/tess_false_alarms/keplerTargets/nearby_stars/"])
    star_file = "/Users/smullally/Science/tess_false_alarms/keplerTargets/target_selection/tic_kic_stars_15.csv"
    star_data = p.read_csv(star_file)
    stars = star_data['ID']
    #Pick up from 100000 and redo 60000+10855, 10882, 37572
    w = 100000+np.array([11162, 11891, 12466, 12470, 14213, 17746])
    complete = parmap.parmap(nearby_kic_make_local, outdir, stars[w], engine='threads')
    
    return complete
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

def get_true_nature(atce, time_range, time_error):
    """
    unkown means does not match anything in the varible files provided.
    """
    
    known_variables_file = "/Users/smullally/Science/tess_false_alarms/keplerTargets/catalogs/varTable.csv"
    vartab = ascii.read(known_variables_file)
    
    #nearby_stars = nearby_kic(atce['star'])
    print(atce['star'][4:])
    nearby_stars = nearby_kic_local(int(atce['star'][4:]), 
                                 "/Users/smullally/Science/tess_false_alarms/keplerTargets/nearby_stars/")
    nearstars = nearby_stars[~nearby_stars['KIC'].mask]
    closekic = nearstars[0]['KIC']
    
    #Get all TCEs associated with the nearby stars    
    nearby_tces = get_tces(nearstars['KIC'], vartab)
    
    status = ""
    
    for neartce in nearby_tces:
        
        periods = [atce['period'], neartce['period']]
        epochs = [atce['epoch'] - atce['epoch_offset'], 
                  neartce['epoch'] - neartce['epoch_offset']]
        
        #time_range = [1711*u.d - exo_const.btjd, 1763*u.d - exo_const.btjd]
        #time_error = atce['duration']*0.5
        
        num, small_diffs = transit_time_match(periods, epochs, time_range, time_error)

        if num>=1:
            if closekic == neartce['kicid']:
                status = ';'.join([status, "PC %s num:%u" % (neartce['signal_name'], num)])
                return status
            else:
                status = ';'.join([status, "CO %s num:%u" % (neartce['signal_name'], num)])

    #import pdb;pdb.set_trace()
    if status == "":
        status = "unknown"
    
    return status
    

#%% 
    
from exovetter import utils
def example():
    #'KIC 10797460'
    time_range = [1711*u.d - exo_const.btjd, 1763*u.d - exo_const.btjd]
        
    atoi = utils.get_mast_tce('K00752.01')
    koi752 = atoi[0]
    koi752['star'] = 'TIC 26539443'
    time_error = koi752['duration']*0.4
    print(koi752)
    status = get_true_nature(koi752, time_range, time_error)
    
    print(status)
    
    atoi = utils.get_mast_tce('K00759.01')
    koi759 = atoi[0]
    koi759['star'] = 'TIC 377782363'
    print(koi759)
    time_error = koi759['duration']*0.5
    status = get_true_nature(koi759, time_range, time_error)
    print(status) #PC
    
    atoi = utils.get_mast_tce('K00759.01')
    koi759 = atoi[0]
    koi759['star'] = 'TIC 377782366'  #fake for a CO test
    time_error = koi759['duration']*0.5
    print(koi759)
    status = get_true_nature(koi759, time_range, time_error)
    print(status)  #CO
    
    atoi = utils.get_mast_tce('K00759.01')
    koi759d = atoi[0]
    koi759d['period'] = 3.3 * u.d
    koi759d['epoch'] = 55024.1 * u.d
    koi759d['star'] = 'TIC 377782366'  #fake for a CO test
    time_error = koi759['duration']*0.5
    print(koi759d)
    status = get_true_nature(koi759d, time_range, time_error)
    print(status) #should be unknown
    
    
#%%
def get_time_range(sector):
    d = dict()
    d['14'] = (1683.3*u.d - exo_const.btjd, 1710.2*u.d - exo_const.btjd)
    d['15'] = (1711.35*u.d - exo_const.btjd, 1737.41*u.d - exo_const.btjd)
    d['25'] = (1983.63*u.d - exo_const.btjd, 2009.30*u.d - exo_const.btjd)
    d['26'] = (2010.26*u.d - exo_const.btjd, 2035.13*u.d - exo_const.btjd)
    
    return d[str(sector)]


import pandas as p
def ephem_match_list():
    """
    Given a file names of ephemerides, create a tce and determine their
    true nature returning a PC, CO or unknown.
    """
    filename = "/Users/smullally/Science/tess_false_alarms/vet_results/joe_mar_2021/first-dataset/qlp_noebplanets_tcesum.csv"    
    #filename = "/Users/smullally/Science/tess_false_alarms/vet_results/joe_mar_2021/second-dataset/spoc_ebplanets_tcesum.csv"

    col_names = ["tic","pn","sector","period","epoch",
                 "depth","dur","snr","disp","reason"]
    
    tce_df = p.read_csv(filename, names=col_names)
    
    outname = filename + '.em.csv'
    dh = open(outname, 'w+')  #switch to 'a' if don't want to overwrite
    
    header = True  #Only writes header on the first loop
    for i, data in tce_df[:].iterrows():
        
        mytce = Tce(period = data['period'] * u.day,
                    epoch = data['epoch'] * u.day,
                    epoch_offset = exo_const.btjd,
                    depth = data['depth'] * exo_const.ppm,
                    duration = data['dur'] * u.hour,
                    snr = data['snr'],
                    disp = data['disp'],
                    ephem_match = "",  #Status to be filled in
                    star = data['tic'],
                    sector = data['sector'],
                    signal_name = "%s-p%u-s%u" % (data['tic'],data['pn'],data['sector']))
        #print(mytce['signal_name'])
        
        time_range = get_time_range(mytce['sector'])
        time_error = mytce['duration'] * 0.4
        
        ephem_match = get_true_nature(mytce, time_range, time_error)
        tce_df.at[i,'ephem_match'] = ephem_match
        #data['ephem_match'] = ephem_match
    

        tce_str = tce_df.iloc[i:i+1].to_csv(None,sep=',', header=header)
        header = False

        dh.write(tce_str)
    
    #outname = filename + '.em.csv'

    #tce_df.to_csv(outname)  
    dh.close()
    return(tce_df)
