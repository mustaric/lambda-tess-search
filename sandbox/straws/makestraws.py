# -*- coding: utf-8 -*-
# Copyright 2017-2018 Orbital Insight Inc., all rights reserved.
# Contains confidential and trade secret information.
# Government Users: Commercial Computer Software - Use governed by
# terms of Orbital Insight commercial license agreement.

"""
Created on Fri Aug 16 19:25:37 2019

TODO
x Save metadata to output directory
x Function to read metadata
o Load straws should figure out if its searching on s3
o Straw names should include sector
o Docstrings
o Investigate async to speed writing the straws
o Add more info to metadata file
o Straw maker may not deal with edges of the ccd correctly
o Is npy the best format for writing straws?
o Do I need to save a time file?

@author: fergal


Concepts
-----------
FFI
    Full Frame Image. A single image from a single camera/ccd taken by
    TESS. Typically, we only want a small postage stamp from each FFI,
    but we want that postage stamp for many FFIs

Datacube
    A 3d numpy array, where each row is an image.

Straw
    A datacube consisting of a small postage stamp image from each an every
    FFI in a sector.

Camera, CCD, col, row
    This 4-tuple uniquely identifies a pixel on the TESS focal plane.

This module contains a class to collect all the imagery from a set of FFIs,
and create straws that tile the entire viewing area, then another class
to construct a datacube from those tiles.

This serves the use case where you want to get all the imagery across all
FFIs for a single star. By dividing up the data this way we can minimise
the time spent on reading files, and downloading data. As such it is an
excellent approach for accessing data on single stars from s3 using an
AWS Lambda.

"""

from __future__ import print_function
from __future__ import division

from pdb import set_trace as debug
import numpy as np
import os

import astropy.io.fits as pyfits
from glob import glob
import json

from common import  METADATA_FILE
from common import makeStrawName


class MakeTessStraw(object):
    def __init__(self, ffiPath, outPath, sector, camera, ccd):
        """
        Inputs
        ----------
        ffiPath
            (str) Path to FFI files on local disk. This path can contain
            FFIs from mulitple ccds or cameras, but can only contain FFIs
            from a single sector.

        outPath
            (str) Location on local disk to store straws
        sector
            (int) sector number to process
        """

        self.outPath = outPath
        self.ffiPath = ffiPath
        self.sector = sector
        self.camera = camera
        self.ccd = ccd
        self.strawSize = 50
        self.midtimes_tbjd = None #Will be filled in later

        #The sector version string is part of the FFI filename
        sectorVersion= {1:120, 3:123}
        try:
            self.dataVersion = sectorVersion[sector]
        except IndexError:
            raise IndexError("sectorVersion string not hardcoded for this sector yet")

        #These must be set in this order
        self.datestampList = self.loadDatestamps()
        self.nColsRows = self.getFfiShape()

        self.do(camera, ccd)

    def loadDatestamps(self):
        """Load a list of datestamps from all the FFIs in `ffiPath`

        Look in the directory `self.ffiPath`, find all the FFI files,
        and extract their datestamps. It would be nice is those were also
        the observation times, but this is the SPOC, so they aren't.

        The datestamps are used for loading the right file from disk
        """
        pattern = "tess*ffic.fits"
        pattern = os.path.join(self.ffiPath, pattern)
        fileList = glob(pattern)
        assert len(fileList) > 0

        f = lambda x: os.path.split(x)[-1].split('-')[0][4:]
        datestamps = map( f, fileList)

        datestamps = sorted(list(set(datestamps)))
        assert len(datestamps) > 0
        return datestamps


    def do(self, camera, ccd):
        nCols, nRows= self.nColsRows

        for i in range(0, nCols, self.strawSize):
            print("Processing column %i" %(i))
            for j in range(0, nRows, self.strawSize):
                straw, times = self.makeStraw(camera, ccd, i, j)
                self.writeStraw(straw, camera, ccd, i, j)
            break

        self.midtimes_tbjd = times
        self.saveMetadata()


    def getFfiShape(self):
        """Get the num cols and rows from the FFI header"""
        ffiName = self.getFfiName(0, self.camera, self.ccd)
        hdr = pyfits.getheader(ffiName, 1)
        nCols = hdr['NAXIS1']
        nRows = hdr['NAXIS2']

        return nCols, nRows

    def makeStraw(self, camera, ccd, col, row):
        """Make a straw at the requested location

        Inputs
        -------------
        camera, ccd, col, row
            (int) Properties of the straw. col and row refer to coordinates of
            the bottom-left corner of the straw.

        Returns
        ---------
        np 3d array
        """
        nCol, nRow = self.strawSize, self.strawSize
        nCadence = len(self.datestampList)
        straw = np.empty( (nCadence, nRow, nCol) )
        midtimes_tbjd = np.empty(nCadence)

        for i in range(nCadence):
            ffiName = self.getFfiName(i, camera, ccd)
            frame, time = self.readFfiSection(ffiName, col, row)
            nr, nc = frame.shape
            straw[i,:nr,:nc] = frame
            midtimes_tbjd[i] = time
        return straw, midtimes_tbjd

    def writeStraw(self, straw, camera, ccd, col, row):
        """
        Write a straw to disk

        Inputs
        -----------
        straw
            (3d npy array) The data to save

        camera, ccd, col, row
            (int) Properties of the straw. col and row refer to coordinates of
            the bottom-left corner of the straw.
        """
        path, fn = makeStrawName(self.outPath,
                                 self.sector,
                                 camera,
                                 ccd,
                                 col,
                                 row)

        if not os.path.exists(path):
            os.makedirs(path)

        fn = os.path.join(path, fn)
        np.save(fn, straw)

    def readFfiSection(self, ffiName, col, row):
        """Read a postage stamp from an FFI file

        This is by far the most expensive function in the class
        and optimisation efforts should focus here
        """

        hdulist = pyfits.open(ffiName, memmap=True)

        img = hdulist[1]
        slCol = slice(col, col + self.strawSize)
        slRow = slice(row, row + self.strawSize)
        data = img.section[slRow, slCol]

        tstart = img.header['TSTART']
        tend = img.header['TSTOP']
        midtime_tbjd = .5 * (tstart + tend)

        hdulist.close()
        return data, midtime_tbjd

    def getFfiName(self, cadenceNum, camera, ccd):
        """Construct the path to an FFI on local disk

        Raises an index error if `cadenceNum` is out of bounds on the
        list of FFIs available
        """

        datestamp = self.datestampList[cadenceNum]

        fn = "tess%s-s%04i-%i-%i-%04i-s_ffic.fits" \
            %(datestamp, self.sector, camera, ccd, self.dataVersion)
        return os.path.join(self.ffiPath, fn)

    def saveMetadata(self):
        """Save a metadata file to a local filestore
        """
        #Convert to a JSON serialisable list
        self.midtimes_tbjd = list(self.midtimes_tbjd)

        fn = os.path.join(self.outPath, METADATA_FILE)
        text = json.dumps(self.__dict__, indent=2)
        with open(fn, 'w') as fp:
            fp.write(text)

