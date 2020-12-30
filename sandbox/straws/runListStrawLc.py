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
#%%
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
    if radius > 3.5:
        radius = 3.5
        
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
for i,ticid in enumerate(tics[10:14]):
    
    jsonInput = {"ticid": str(int(ticid)), 
                 "straw_bucket": "tess-straws", 
                 "lc_bucket": "straw-lightcurves", 
                 "ap_radius": str(radii[i+10]), 
                 "sector": "16" }
    print(jsonInput)
    client.invoke(FunctionName="strawLcMaker",
                  InvocationType="Event",
                  Payload = json.dumps(jsonInput)
                  )
b=time.time()

total=b-a
print(total)
#%%
#This list is for a larger S0001 run.
ticfile = "/Users/smullally/TESS/lambdaSearch/s01StrawRun/s0001-cam1-Tmag11-14.csv"

ticlist = p.read_csv(ticfile, delimiter=',', comment='#')
tics= ticlist['ID']
Tmags = ticlist['Tmag']
radii= np.array(list( map (lambda mag:  pick_radius(mag), Tmags)))
#%%
client = boto3.client('lambda')

a= time.time()
for i in range(4000,6000):
    print("TIC: %u  Tmag: %f   radii: %f " % (tics[i], Tmags[i], radii[i]))
    jsonInput = {"ticid": str(int(tics[i])), 
                 "straw_bucket": "tess-straws", 
                 "lc_bucket": "straw-lc-timing-test",
                 #"lc_bucket": "straw-lightcurves", 
                 "ap_radius": str(radii[i]), 
                 "sector": "1" }
    #print(jsonInput)
    client.invoke(FunctionName="strawLcMaker",
                  InvocationType="Event",
                  Payload = json.dumps(jsonInput)
                  )
b=time.time()

total=b-a
print(total)