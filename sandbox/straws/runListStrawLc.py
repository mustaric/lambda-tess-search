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

ticfile = "/Users/smullally/TESS/lambdaSearch/strawTests/s0001-1-1_targlist1.txt"

ticlist = np.loadtxt(ticfile)

client = boto3.client('lambda')

a= time.time()
for tic in ticlist[80:101]:
    
    jsonInput = {"ticid": str(int(tic)), 
                 "straw_bucket": "tess-straws", 
                 "lc_bucket": "straw-lightcurves", 
                 "ap_radius": "2", 
                 "sector": "1" }
    print(jsonInput)
    client.invoke(FunctionName="strawLcMaker",
                  InvocationType="Event",
                  Payload = json.dumps(jsonInput)
                  )
b=time.time()

total=b-a
print(total)