#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is code to create a list of targets observed by both 
Kepler original mission and TESS Sectors 14/15
We want ones that don't have known EBs or large planets.
So removing those in planet and EB catalogs.

Created on Tue Jun 30 15:19:47 2020

@author: smullally
"""


#Read in different catalogs into Pandas dataframes.
#Implement selection criteria.
#do not want targets with known EBs or exoplanets.
#At the very least, put those in a separate list.

#Want targets brighter than 14.

import pandas as p
from astropy.io import ascii
from astroquery.mast import Catalogs
import numpy as np
from astropy.io import ascii


#%%
ddir = "/Users/smullally/Science/tess_false_alarms/keplerTargets/catalogs/"

#List of Kepler IDs that we know something about.
keplerstarfile = ddir + "kepler_steller_17.csv"
kepstar = p.read_csv(keplerstarfile,sep="|")
bright_kepstar = kepstar['kepmag'] < 14


#planets
koifile = ddir + "cumulative.csv"
data = ascii.read(koifile)
kois = data.to_pandas()

#EBs
ebfile = ddir + "villanova_eb_catalog.csv"
ebs = p.read_csv(ebfile,sep=',',comment="#")

#%%
# Now I need kicids that are not in the other lists.

clean_stars = set(kepstar[bright_kepstar]['kepid']) - set(kois['kepid'])
clean_stars = set(clean_stars) - set(ebs['KIC'])
print("Number of clean_stars %u" % (len(clean_stars)))

#%%
#Get TIC ID information for good_stars

stars = np.array(list(clean_stars))

#tic_data = Catalogs.query_criteria(catalog="Tic", KIC=stars[0:20000])
#ascii.write(tic_data, '/Users/smullally/Python_Code/lambda-tess-search/falseAlarms/target_selection/tic_kic_stars1.csv',format='csv')

#tic_data = Catalogs.query_criteria(catalog="Tic", KIC=stars[20000:40000])
#ascii.write(tic_data, '/Users/smullally/Python_Code/lambda-tess-search/falseAlarms/target_selection/tic_kic_stars2.csv',format='csv')

#tic_data = Catalogs.query_criteria(catalog="Tic", KIC=stars[40000:])
#ascii.write(tic_data, '/Users/smullally/Python_Code/lambda-tess-search/falseAlarms/target_selection/tic_kic_stars3.csv',format='csv')

#This can take a while
tic_data = Catalogs.query_criteria(catalog="Tic", KIC=stars)
filename = '/Users/smullally/Science/tess_false_alarms/keplerTargets/target_selection/tic_kic_stars.csv'
ascii.write(tic_data, filename ,format='csv')

#%%

targets = tic_data
#Creat code to read in the csv file. TODO
targets = ascii.read(filename, format='csv')

#%%
#Depends on radial distance (2 pixels?) and magnitude differrence (within 2 magnitudes?)
#
def aperture_radii(Tmag):
    
    #See Stassun 2018, AJ 156 Section 3.3.3
    if Tmag < 4:
        Tmag = 4
    Npix = -0.2592 * Tmag**3 + 7.7410*Tmag**2 - 77.7918*Tmag + 274.2898
    
    if Npix < 0:
        Npix = 2
    
    return Npix * 20.25
#%%
#This take over an hour ot run.
dist = 20.25*10 # radial Distance in arcseconds 20.25 is size of pixel.
dmag = 2 # count number of stars brighter than Tmag-dmag of your star.

for star in targets[~np.isnan(targets['contratio'])][20:]:
    
    starid = star['ID']
    name = "TIC %s" % starid
    To = star['Tmag']
    Ao = aperture_radii(To)/2
    
    print(star['ID', 'Tmag'], Ao)
    near_stars = Catalogs.query_object(objectname=name, radius=dist/3600, catalog="Tic")
    bright = near_stars['Tmag']< To + dmag
    
    #print(near_stars[bright]['ID','Tmag','dstArcSec'])
    
    count = 0
    
    for near in near_stars[bright][1:]:
        Astar = aperture_radii(near['Tmag'])/2
        #print(Astar)
        
        if (Astar + Ao) > near['dstArcSec']:
            count += 1
            
    print(count)
    if count == 0:
        targets['contratio'][targets['ID'] == starid] =  -1  
#%%
newfile = "/Users/smullally/Science/tess_false_alarms/keplerTargets/target_selection/tic_kic_stars_contr.csv"
ascii.write(targets, newfile,format='csv')
    
#%%
#Find uncrowded stars
contamination_limit = 0.15
uncrowded = targets['contratio'] < contamination_limit

good_targets = targets[uncrowded]
print(len(good_targets))

#%%
#2020-10-09
#Use tess-point to determine which stars are actually on silicon.
from tess_stars2px import tess_stars2px_function_entry
import pandas as p

ticids = list(good_targets['ID'])
ras = list(good_targets['ra'])
decs = list(good_targets['dec'])

# do one first
outid, outlong, outlat, sector, cam, ccd, colpix, rowpix, scinfo \
    =  tess_stars2px_function_entry(ticids[0], ras[0], decs[0])
#%%
# Now do all in a list, sending in the scinfo from above for speed.
outid, outlong, outlat, sector, cam, ccd, colpix, rowpix, scinfo \
   = tess_stars2px_function_entry(ticids[0:2000], ras[0:2000], decs[0:2000])

data = {'TIC':outid, 'Sector':sector, 'Camera':cam, 'ccd':ccd, \
        'colpix':colpix, 'rowpix':rowpix} 

obs = p.DataFrame(data = data)

obsfile = "/Users/smullally/Science/tess_false_alarms/keplerTargets/target_selection/tic_uncrowded_bright-sectorsObs.csv"

obs.to_csv(obsfile)
    
#good_targets['sector'] = sector

