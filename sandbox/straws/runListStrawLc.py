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
for i,ticid in enumerate(tics):
    
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