# -*- coding: utf-8 -*-
# Copyright 2017-2018 Orbital Insight Inc., all rights reserved.
# Contains confidential and trade secret information.
# Government Users: Commercial Computer Software - Use governed by
# terms of Orbital Insight commercial license agreement.

"""
Created on Tue Oct 22 21:22:36 2019

@author: fergal
"""

from __future__ import print_function
from __future__ import division

import boto3

import numpy as np
import json
import os
import io

import common

class LoadTessCube(object):
    """
    Load a datacube of TESS imagery from straws stored on disk.
    """

    def __init__(self, path, sector, camera, ccd, checkver=True):
        """

        Parameters
        ----------
        path : TYPE
            DESCRIPTION.
        sector : TYPE
            DESCRIPTION.
        camera : TYPE
            DESCRIPTION.
        ccd : TYPE
            DESCRIPTION.
        checkver : bool, optional
            Setting checkver to False turns off the check that the 
            version of the metadata file agrees with the version
            of the code. This should only be done for debugging purposes.


        Returns
        -------
        None.

        """

        #Set path to None for some testing
        if path is not None:
            self.path = path
            self.sector = sector
            self.camera = camera
            self.ccd = ccd
            self.loadMetadata()

    def __repr__(self):
        return "<TessCube object for sector %i, camera %i, ccd %i. Data at %s>" \
            %(self.sector, self.camera, self.ccd, self.path)

    def __call__(self, col, row):
        return self.get(col, row, min_size_pix=20)

    def loadMetadata(self):
        """Load metadata on the straws stored in `path`

        Metadata is stored in a json file and contains details like ccd sizes,
        number of cadences, strawsize, etc.
        """
        fn = common.getMetadataPath(self.path, 
                                    self.sector,
                                    self.camera,
                                    self.ccd)

        with open(fn) as fp:
            props = json.load(fp)
            
        self.checkMetadataSanity(props)
        self.setMetadataFromDict(props)
        
    def setMetadataFromDict(self, props):
        self.__dict__.update(props)
        self.nCols, self.nRows = self.nColsRows
        self.nCadences = len(self.datestampList)

    def checkMetadataSanity(self, props):
        try:
            dataver = props['__straw_version__']
        except KeyError:
            raise ValueError("Old obsolete metadata file without version info found")
            
        expectver = common.STRAW_VERSION
        if  dataver !=  expectver:
            raise ValueError("Expected version %s straws, got version %s" %(expectver, dataver))

        assert self.sector == props['sector']
        assert self.camera == props['camera']
        assert self.ccd == props['ccd']
        
    def getMidTimestamps(self):
        """Return the cadence mid times as stored in the metadata
        
        See make straws for the details of how this value is calculated
        """
        
        try:
            timestamps = self.midtimes_tbjd
        except AttributeError:
            raise AttributeError("metadata doesn't contain timestamps")
            
        return np.array(timestamps)

    def getRelativeCadenceNumbers(self):
        """Return a integers from zero to length of datacube"""
        return np.arange(self.nCadences, dtype=int)
    
    def get(self, col, row, min_size_pix=None):
        """Get a data cube

        The data cube is garaunteed to be square and at least `min_size_pix`
        on a side. However, because it constructs that datacube whose bounding
        box aligns with the straws its reading data from, the actual size
        may be larger than `min_size_pix`, and the requested (`col`, `row`)
        may not be at the centre of the image.

        Inputs
        -------------
        col, row
            (int) Location on CCD to load a straw for

        Optional Inputs
        -----------------
        min_size_pix
            (int) Minimum width and height of the returned datacube

        Returns
        -----------
        cube
            (np 3d array) of shape (nCadence, nRows, nCols)
        target_col, target_row
            (float) The index in `image` corresponding to (`col`, `row`).
            For example, if the request is for a 30x30 pixel stamp around
            the postion cr= 301, 602, the resulting target_col, _row might be
            (1,2)
        """
        if min_size_pix is None:
            min_size_pix = self.strawSize

        c0, c1, r0, r1 = self.pickBbox(col, row, min_size_pix)
        colSize = c1 - c0
        rowSize = r1 - r0
        image = np.empty( (self.nCadences, rowSize, colSize) )

        ds = self.strawSize
        for i in range(c0, c1, ds):
            for j in range(r0, r1, ds):
                straw = self.getStraw(i, j)
                
                assert straw.shape == (self.nCadences, ds, ds)

                dCol = i - c0
                dRow = j - r0
                sc = slice(dCol, dCol + ds)
                sr = slice(dRow, dRow + ds)
                image[:, sr, sc] = straw

        target_col = col - c0
        target_row = row - r0
        return image, target_col, target_row

    def pickBbox(self, col, row, size_pix):
        """Pick the bounding box around (col, row) for the returned data cube

        The bounding box will be

        * square
        * The width will be > `size_pix`
        * The width will be an integer times the `strawSize`

        Inputs
        -------
        col, row
            (float) Location of centre of region of interest
        size_pix
            (int) Minimum size of returned bounding box. The bounding box
            will probably be bigger than this request.


        Returns
        ----------
        4-tuple of col and row values defining the bounding box.
        """
        if not self.isInBounds(col, row):
            raise ValueError("Requested col,row (%g, %g) is out of bounds" %(col, row))
        assert(size_pix > 0)

        ds = .5 * size_pix
        c0 = common.roundToNearestBelow(max(col-ds, 0), self.strawSize)
        c1 = common.roundToNearestAbove(min(col+ds, self.nCols), self.strawSize)

        r0 = common.roundToNearestBelow(max(row-ds, 0), self.strawSize)
        r1 = common.roundToNearestAbove(min(row+ds, self.nRows), self.strawSize)

        return c0, c1, r0, r1

    def isInBounds(self, col, row):
        """Test if the requested col,row actually fall on disk

        Inputs
        -------------
        col, row
            (int)

        Returns
        ----------
        boolean
        """

        if col < 0 or col >= self.nCols:
            return False

        if row < 0 or row >= self.nRows:
            return False

        return True

    def getStraw(self, col, row):
        """ Load a straw from disk

        Inputs
        -------------
        camera, ccd, col, row
            (int) Properties of the straw. col and row refer to coordinates of
            the bottom-left corner of the straw.

        """
        longPath, fn = common.makeStrawName(self.path,
                                 self.sector,
                                 self.camera,
                                 self.ccd,
                                 col,
                                 row)

        straw = self.loadStrawFromUri(longPath, fn)
        return straw

    def loadStrawFromUri(self, strawPath, fn):
        if not os.path.exists(strawPath):
            raise IOError("Path %s not found" %(strawPath))

        fn = os.path.join(strawPath, fn)
        if not os.path.exists(fn):
            raise IOError("File %s not found" %(fn))

        return np.load(fn)



class LoadTessCubeS3(LoadTessCube):
    """Load straws from S3 instead of a local disk"""

    def __init__(self, bucket, path, sector, camera, ccd, region='us-east-1'):
        #bucket is a string. self.bucket is an object
        self.bucketName = bucket
        self.s3 = boto3.resource('s3', region_name=region) 
        self.path = path
        self.sector = sector
        self.camera = camera
        self.ccd = ccd
        self.loadMetadata()

    def loadStrawFromUri(self, strawPath, fn):
        #boto stuff goes here
        uri = os.path.join(strawPath, fn)
        obj = self.s3.Object(self.bucketName, uri)
        thebytes = obj.get()['Body'].read()
        return np.load(io.BytesIO(thebytes))

    def loadMetadata(self):
        """Load metadata on the straws stored in `path`

        Metadata is stored in a json file and contains details like ccd sizes,
        number of cadences, strawsize, etc.
        """
        uri = common.getMetadataPath(self.path, 
                                    self.sector,
                                    self.camera,
                                    self.ccd)
        print(uri)
        obj = self.s3.Object(self.bucketName, uri)
        print(obj)
        thebytes = obj.get()['Body'].read()

        props = json.loads(thebytes)
        self.checkMetadataSanity(props)
        self.setMetadataFromDict(props)
