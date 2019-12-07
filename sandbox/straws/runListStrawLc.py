#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 21:08:24 2019

@author: smullally
"""

import boto3
import numpy as np
import json
import time
import pandas as p

ticfile = "/Users/smullally/TESS/lambdaSearch/strawTests/kdwarf_cam1_input.csv"

ticlist = p.read_csv(ticfile)
tics= ticlist['ticid']
radii= ticlist['radii']
#%%
client = boto3.client('lambda')

a= time.time()
for i,ticid in enumerate(tics[1:10]):
    
    jsonInput = {"ticid": str(int(ticid)), 
                 "straw_bucket": "tess-straws", 
                 "lc_bucket": "straw-lightcurves", 
                 "ap_radius": str(radii[i]), 
                 "sector": "1" }
    print(jsonInput)
    client.invoke(FunctionName="strawLcMaker",
                  InvocationType="Event",
                  Payload = json.dumps(jsonInput)
                  )
b=time.time()

total=b-a
print(total)


#%%
def pick_radius(Tmag):
    radius =(-0.34 * Tmag) + 5.6
    if radius < 0.80:   #Smallest aperture is on average 2 pixels ~ 14th magnitude
        radius = 0.80
    if radius > 3.0:
        radius = 3.0
        
    return radius
    
        
#%%
#This list is from sector 16
ticfile = "/Users/smullally/TESS/lambdaSearch/s16StrawRun/s0016_cam3_dwarfStars.csv"

ticlist = p.read_csv(ticfile, delimiter=',', comment='#')
tics= ticlist['ID']
Tmags = ticlist['Tmag']
radii= np.array(list( map (lambda mag:  pick_radius(mag), Tmags)))
#%%
client = boto3.client('lambda')

a= time.time()
for i,ticid in enumerate(tics[0:2]):
    
    jsonInput = {"ticid": str(int(ticid)), 
                 "straw_bucket": "tess-straws", 
                 "lc_bucket": "straw-lightcurves", 
                 "ap_radius": str(radii[i]), 
                 "sector": "16" }
    print(jsonInput)
    client.invoke(FunctionName="strawLcMaker",
                  InvocationType="Event",
                  Payload = json.dumps(jsonInput)
                  )
b=time.time()

total=b-a
print(total)