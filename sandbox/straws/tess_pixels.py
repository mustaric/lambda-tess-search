#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 12:13:33 2019

@author: smullally
"""

from astropy.wcs import WCS
from astroquery.mast import Observations
from astroquery.mast import Catalogs
from astropy.coordinates import SkyCoord
from astroquery.mast import Tesscut
from astropy.io import fits
import boto3
import io

def id_wcs_file(secid):
   """
   Use astroquery to locate the FFI file that contains the WCS.
   Download it
   Return the filtered data product.
   """
    
   obs_table = Observations.query_criteria(obs_id=secid)
   products = Observations.get_product_list(obs_table)
   filtered = Observations.filter_products(
               products, productSubGroupDescription="FFIC",
               mrp_only=False)

   return filtered
    

def download_wcs_file(filtered, local_dir, n=10, cloud=False):
    """
    filtered is results if id_wcs_file
    local_dir is the location to put the file
    cloud determines whether the file should first be pulled from the cloud.
    
    For lambda local_dir should be /tmp
    """
    
    if cloud:
        Observations.enable_cloud_dataset(provider='AWS')
            
    obsslice = slice(n,n+1)
    manifest=Observations.download_products(filtered[obsslice], 
                                            download_dir=local_dir, 
                                            mrp_only=False)
    
    return manifest['Local Path'][0]
     

def make_secid(sector, camera, ccd):

    secid = "tess-s%04u-%u-%u" % (sector, camera, ccd)   
    
    return secid    
        
def get_xy(filename, coord):
    """
    Input:
        local file name of the FFI containing the WCS.
        coord of the star
    return:
        xpos: column position
        ypos: row position
    """
    
    hdr = fits.getheader(filename, ext=1)
    w = WCS(hdr)
    pix = w.all_world2pix(coord.ra.deg, coord.dec.deg, 0)
    xpos = float(pix[0])
    ypos = float(pix[1])
    print(xpos,ypos)
    
    return(xpos, ypos)


def ticid_to_coord(ticid):
    """
    Input:
        ticid: string containing the TIC Identificaiton number.
    Returns:
        coord: astropy SkyCoordinate object of ra and dec
    """
    info = Catalogs.query_criteria(catalog='Tic', ID=ticid)
    ra = info['ra'][0]  # deg
    dec = info['dec'][0]  # deg
    coord = SkyCoord(ra, dec, unit='deg')   
        
    return coord
        
def get_object_coords(ticid, sector, nFFI=10, cloud = False, local_dir = "."):
    """
    Retrieve the camera, ccd, column row of the ticid of interest.
    Downloads a single FFI image in order to retrieve the FITS header.
    Input:
        ticid: tess ID
        sector : integer sector number
        nFFI: indicates which FFI cadence to use.
        cloud: False. Set to true if running on AWS
        local_dir: Is the location to store the FFI image before reading header.
        
    Returns:
        camera, ccd, column, row
        
        if they are all zero, then there is no data to be had
    """
    
    coord = ticid_to_coord(ticid)
    results = Tesscut.get_sectors(coord,radius=0)
    
    camera = 0
    ccd = 0
    col = 0
    row = 0
    
    if len(results)>0:
        
        want = results['sector'] == sector
        if len(results[want]) > 0:
            
            camera = results[want]['camera'][0]
            ccd = results[want]['ccd'][0]

    if camera > 0:
        secid = make_secid(sector, camera, ccd)
        ffi_file_path = get_wcsfile(secid)
        #filtered = id_wcs_file(secid)

    #ffi_file_path = download_wcs_file(filtered, local_dir, n=nFFI, cloud=cloud)
    #I need to create a file lookup in here, but for the moment.

    col, row = get_xy(ffi_file_path, coord)

            
    return camera, ccd, col, row

def get_wcsfile(secid):
    """
    return wcs header
    """
    bucket = "tess-straws"
      
    ffi_file_name = "ffiwcs/%s_wcs.fits" % secid[-9:]

    print(secid)
    s3 = boto3.resource('s3')
    obj = s3.Object(bucket, ffi_file_name)
    thebytes = obj.get()['Body'].read()
    
    return io.BytesIO(thebytes)
    