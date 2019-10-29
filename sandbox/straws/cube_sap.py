#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 10:29:58 2019

@author: smullally
"""

import extractlc as exlc
import numpy as np

def get_fluxes(cube, centroid=None, radius_pix=3):
    """
       Input:
       cube: 3d array of flux images with dimensions (time, column, row)
       
       centroid is (col, row) of the location of the star in the cube.
       
       radius_pix is the size of the circular aperture to use in pixels.
       
       Note this code uses the entire cube to deterine the sky. 
       The aperture photometry is done using a circle of the specified radius.
       The center of the pixel needs to be <= the distance to be considered
       in the aperture.
       
       Returns:
       fluxMinusSky: Simple, aperture photometry flux minus the measured sky flux
       sky: the measured sky
       Both are in the same units as the individual pixels in the cube.
    """

    sky = exlc.measureEmptySky(cube)
    
    #Use center if centroid is not specified
    if centroid is None:
        centroid = np.array(cube.shape[1:]) / 2.

    #replace nans with 0.0
    if ~np.all(np.isfinite(cube.sum(axis=0))):
        cube = np.nan_to_num(cube)

    #Choose a simple circular aperture for the star.  
    row,col = np.shape(cube)[1:3]  
    aper = exlc.chooseSimpleAperture((col,row), starloc=centroid, 
                                     radius=radius_pix)

    fluxMinusSky = exlc.performAperturePhotometry(cube, sky, aper)

    ave_image = np.nanmean(cube, axis=0)
    
    return fluxMinusSky, sky, ave_image