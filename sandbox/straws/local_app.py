#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 21:05:24 2019

@author: smullally
"""

#Main file that runs locally on my machine for developing 
#the straw light curve generation.


import loadstraws as ls
import numpy as np
import tessid_to_pix as tesspx
import pixels_to_lightcurve as ptl
import write_lightcurve


def lambda_handler(event, context):
    
    ticid = event['ticid']
    straw_bucket = event['straw_bucket']
    lc_bucket = event['lightcurves']
    ap_radii = event['ap_radii']
    
    #length is the number of unique sectors/cam/ccd for ticid
    sectors, camera, ccds, ra, dec = tesspx.get_pixels(ticid)
    
    for s in sectors:
        
        #This is written
        path = straw_bucket + "/s%0u" % (s)
        cubeObj = ls.LoadTessCube(path)
        cube, cube_col, cube_row = cubeObj.get(camera, ccd, col, row, min_size_pix = 40)
        #midtime = cubeObj.getMidTimestamps()
        midtime = cubeObj.getRelativeCadenceNumbers()
        
        flux, bkg = ptl.do_sap(cube, cube_col, cube_row, ap_radii)
        
        write_lightcurve.to_fits_local(lc_bucket,ticid, sector, camera, ccd)
    
    
    return {
        "statusCode": 200,
        "body": json.dumps({
                "number": str(len(sectors))
                })
        }

def test1():
    #Bucket names should be local file directories for the moment.
    event = {"ticid": "", "straw_bucket": "tess-straws", 
             "lc_bucket":"lightcurves", "ap_radii:":"2.0"}
    context = {}
    val = lambda_handler(event,context)
    print(val)
