#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 16:24:45 2019

@author: smullally
"""

from astropy.table import Table
from datetime import datetime
from astropy.io import fits

def to_fits_local(writepath, output):
    """
    Write basic information to a fits file.
    time, sap_flux, background
    Header contains sector, camera, ccd
    """
    ticid = output['ticid']
    lc_meta = {
            'TELESCOP': 'TESS',
            'CAMERA': output['camera'],
            'SECTOR': output['sector'],
            'CCD': output['ccd'],
            'OBJECT': f'TIC {ticid}',
            'RADESYS': 'ICRS',
            'AP_RAD': output['ap_radius'],
            'CUBESHP': str(output['cube'].shape),
            'CUBEROW': output['cube_row'],
            'CUBECOL': output['cube_col'],
            'ROW': output['row'],
            'COL': output['col'],
            'DATE': datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')}

        # f4 = np.float32, f8 = np.float64, i4 = np.int32
    lc_tab = Table(names=('TIME', 'SAP_FLUX', 'SAP_BKG'),
                   dtype=('f8', 'f4', 'f4'),
                   data=(output['midtime'], output['sap_flux'], output['bkg']),
                   meta=lc_meta)
    c1 = fits.Column(name='TIME', array=output['midtime'], format='f8')
    c2 = fits.Column(name='SAP_FLUX', array = output['sap_flux'], format='f4')
    c3 = fits.Column(name='SAP_BKG', array = output['bkg'], format='f4')
    
    t = fits.BinTableHDU.from_columns([c1, c2, c3], name = 'Light Curve')
    
    
    image = output['av_image']
    
    primary_hdr = fits.PrimaryHDU(header=fits.Header(lc_meta))
    hduim = fits.ImageHDU(image, name = 'Average Image')
    
    hdul = fits.HDUList([primary_hdr, t, hduim])
    
    
    sector = output['sector']
    camera = output['camera']
    ccd = output['ccd']
    sec_id = f's{sector:04}-{camera}-{ccd}'
    basename = f'tic{ticid:0>12}_{sec_id}_lcc.fits'
    outfilename = writepath + basename
    
    #lc_tab.write(outfilename, format='fits')
    hdul.writeto(outfilename, overwrite=True)
    
    
    return outfilename