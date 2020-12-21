#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 09:18:01 2020

@author: smullally

Generate light curves. Each function is a different way to generate a light curve.

"""


import eleanor
import lightkurve as lk


def eleanor_pca(ticid, sector, pxsize = 19):
    
    star = eleanor.Source(tic = ticid, sector = sector)
    
    data = eleanor.TargetData(star, height=pxsize, width=pxsize, bkg_size=31, 
                              do_psf=False, do_pca=True)
    
    return data.time, data.pca_flux, data.quality

#Could do a single sector FFI light curve wiih lightkurve

def eleanor_corr(ticid, sector, pxsize = 19):
    
    star = eleanor.Source(tic = ticid, sector = sector)
    
    data = eleanor.TargetData(star, height=pxsize, width=pxsize, bkg_size=31, 
                              do_psf=False, do_pca=True)
    
    return data.time, data.corr_flux, data.quality


def hlsp(ticid, sector, author="tess-spoc"):
    """
    

    Parameters
    ----------
    ticid : int
        DESCRIPTION.
    sector : int
        Sector of observations to vet
    author : string, OPTIONAL
        options include tess-spoc and tess-qlp.
        The default is "tess-spoc".

    Returns
    -------
    lc : lightkurve object
        lightkurve object of the data requested.

    """
    
    #print(f'TIC {ticid}')
    
    lc = lk.search_lightcurve(f"TIC {ticid}", sector=sector,
                              cadence="ffi",author=author).download()
    
    return lc
    
    