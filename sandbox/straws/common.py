# -*- coding: utf-8 -*-
# Copyright 2017-2018 Orbital Insight Inc., all rights reserved.
# Contains confidential and trade secret information.
# Government Users: Commercial Computer Software - Use governed by
# terms of Orbital Insight commercial license agreement.

"""
Created on Tue Oct 22 21:27:20 2019

@author: fergal
"""

from __future__ import print_function
from __future__ import division

from pdb import set_trace as debug
import numpy as np

import os

METADATA_FILE = 'metadata.json'

def makeStrawName(path, sector, camera, ccd, col, row):
    """Constructh the straw name

    Does no testing to ensure the input parameters are sane

    Inputs
    ---------
    path
        (str) Location of file on disk/s3
    sector, camera, ccd, col, row,
        (int) Parameters of the straw file

    Returns
    -----------
    A string
    """
    path = os.path.join(path,
                        "sector%02i" %(sector),
                        "camera%i" %(camera),
                        "ccd%i" %(ccd))

    fn = "straw-%i-%i-%03i-%03i.npy" %(camera, ccd, col, row)
    return path, fn


def roundToNearestBelow(x, level):
    """Round a number down to the nearest round number

    e.g roundToNearestBelow(180, 100) --> 100
    """

    level = float(level)
    val = np.floor(x / level) * level
    return val.astype(int)

def roundToNearestAbove(x, level):
    """Round a number up to the nearest round number

    e.g roundToNearestBelow(120, 100) --> 200
    """

    level = float(level)
    val = np.ceil(x / level) * level
    return val.astype(int)


def getValuesInRange(x, dx, minx, maxx):
    """Find x1, x2, obeying certain conditions

    The conditions are

    * 0 <= x0 <= x <= x1 <= maxx
    * x0 + dx == x1
    """

    d2 = int(dx / 2.)

    x0 = max(x - d2, minx)
    x1 = min(x + dx, maxx)
    x0 = x1 - dx
    return x0, x1

