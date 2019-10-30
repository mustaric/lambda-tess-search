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
import tess_pixels as tesspx
import cube_sap 
import write_lightcurve
import json
import boto3

def lambda_handler(event, context):
    
    ticid = event['ticid']
    straw_bucket = event['straw_bucket']
    lc_bucket = event['lc_bucket']
    ap_radius = float(event['ap_radius'])
    sector = int(event['sector'])
    
    #Get location of the star for this sector.
    #For cloud use cloud = True and local_dir = /tmp
    camera, ccd, col, row = tesspx.get_object_coords(ticid, sector, \
                                        nFFI=50, cloud = True, local_dir = "/tmp")
     
    
    #Retrieve the Cube.
    cubeObj = ls.LoadTessCubeS3(straw_bucket, "", sector)
    cube, cube_col, cube_row = cubeObj.get(camera, ccd, col, row, 
                                           min_size_pix = 40)

    midtime = cubeObj.getMidTimestamps()
   
    try:
        quality = cubeObj.getQualityFlags()
    except:
        quality = np.zeros(len(midtime))
    sap_flux, bkg, av_image = cube_sap.get_fluxes(cube, centroid=(cube_col, cube_row), 
                               radius_pix=ap_radius)

    
    output  =locals()
    cube_shape=np.shape(cube)
    writepath = "/tmp/"
    filename, basename = write_lightcurve.to_fits_local(writepath, output)
    s3_client = boto3.client('s3')
    esp = s3_client.upload_file(filename, lc_bucket, basename)
    
    
    return {
        "statusCode": 200,
        "body": json.dumps({
                "filename": filename,
                "basename": basename,
                "cubesize": str(cube.shape)
                })
        }

def test1():
    #Bucket names should be local file directories for the moment.
    event = {"ticid": "147424478", 
             "straw_bucket": "/Users/smullally/TESS/tess-straws/", 
             "lc_bucket":"lightcurves", 
             "ap_radius":"3", 
             "sector":"1"}
    context = {}
    val = lambda_handler(event,context)
    print(val)

def test2():
    #osciallatory star
    event = {"ticid": "147424426", 
             "straw_bucket": "tess-straws", 
             "lc_bucket":"straw-lightcurves", 
             "ap_radius":"2", 
             "sector":"1"}
    context = {}
    val = lambda_handler(event,context)
    print(val)
    
def test_loadstraws3():
        
    #Faking it for testing
    sector = 1
    camera = 1
    ccd = 1
    col = 227.5
    row = 255.1
    path = ""
    bucket = "tess-straws"
    cubeObj = ls.LoadTessCubeS3(bucket, path, sector)
    cube, cube_col, cube_row = cubeObj.get(camera, ccd, col, row, 
                                           min_size_pix = 40)
    print(cube.shape)
    