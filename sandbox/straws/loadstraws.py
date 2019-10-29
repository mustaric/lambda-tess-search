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

from pdb import set_trace as debug
from boto.s3.connection import S3Connection

import numpy as np
import json
import os
import io

import common

class LoadTessCube(object):
    """
    Load a datacube of TESS imagery from straws stored on disk.
    """

    def __init__(self, path):

        #Set path to None for some testing
        if path is not None:
            self.path = path
            self.loadMetadata()

    def __repr__(self):
        return "<TessCube object for sector %s. Data at %s>" %(self.sector, self.path)

    def __call__(self, camera, ccd, col, row):
        return self.get(camera, ccd, col, row, 20)

    def loadMetadata(self):
        """Load metadata on the straws stored in `path`

        Metadata is stored in a json file and contains details like ccd sizes,
        number of cadences, strawsize, etc.
        """
        fn = self.getMetadataUrl()
        with open(fn) as fp:
            props = json.load(fp)

        self.__dict__.update(props)
        self.nCols, self.nRows = self.nColsRows
        self.nCadences = len(self.datestampList)

    def getMidTimestamps(self):
        """Return the cadence mid times as stored in the metadata
        
        See make straws for the details of how this value is calculated
        """
        
        try:
            timestamps = self.midtimes_tbjd
        except AttributeError:
            raise AttributeError("metadata doesn't contain timestamps")
            
        return np.array(timestamps)

    def getQualityFlags(self):
        """Return the quality flags from the headers.
        """
        
        try:
            qflags = self.qualityFlags
        except AttributeError:
            raise AttributeError("metadata doesn't contain qualityFlags")
            
        return np.array(qflags)


    def getRelativeCadenceNumbers(self):
        """Return a integers from zero to length of datacube"""
        return np.arange(self.nCadences, dtype=int)
    
    def get(self, camera, ccd, col, row, min_size_pix=None):
        """Get a data cube

        The data cube is garaunteed to be square and at least `min_size_pix`
        on a side. However, because it constructs that datacube whose bounding
        box aligns with the straws its reading data from, the actual size
        may be larger than `min_size_pix`, and the requested (`col`, `row`)
        may not be at the centre of the image.

        Inputs
        -------------
        camera, ccd, col, row
            (int) Properties of the straw. col and row refer to coordinates of
            the bottom-left corner of the straw.

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
                straw = self.getStraw(camera, ccd, i, j)

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

    def getStraw(self, camera, ccd, col, row):
        """ Load a straw from disk

        Inputs
        -------------
        camera, ccd, col, row
            (int) Properties of the straw. col and row refer to coordinates of
            the bottom-left corner of the straw.

        """
        longPath, fn = common.makeStrawName(self.path,
                                 self.sector,
                                 camera,
                                 ccd,
                                 col,
                                 row)



        straw = self.loadStrawFromUri(longPath, fn)
        return straw

    def getMetadataUrl(self):
        return os.path.join(self.path, common.METADATA_FILE)

    def loadStrawFromUri(self, path, fn):
        if not os.path.exists(path):
            raise IOError("Path %s not found" %(path))

        fn = os.path.join(path, fn)
        if not os.path.exists(fn):
            raise IOError("File %s not found" %(fn))

        return np.load(fn)



def LoadTessCubeS3(LoadTessCube):
    """Load straws from S3 instead of a local disk"""

    def __init__(self, bucket, path, aws_key, aws_secret):

        #bucket is a string. self.bucket is an object
        self.bucketName = bucket
        self.aws_secret = aws_secret
        self.aws_key = aws_key

        self.aws_connection = S3Connection(aws_key, aws_secret)
        self.bucket = self.aws_connection.get_bucket(bucket)

        self.loadMetadata(bucket, path)

    def getMetadataUri(self):
       return os.path.join("s3://", self.bucket, self.path, common.METADATA_FILE)

    def loadStrawFromUril(self, path, fn):
        #boto stuff goes here
        uri = os.path.join(path, fn)
        strr = self.bucket.get_key(uri).get_contents_as_string()
        return np.load(io.StringIO(strr))

