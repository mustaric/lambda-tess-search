# -*- coding: utf-8 -*-
# Copyright 2017-2018 Orbital Insight Inc., all rights reserved.
# Contains confidential and trade secret information.
# Government Users: Commercial Computer Software - Use governed by
# terms of Orbital Insight commercial license agreement.

"""
Created on Fri Oct 25 22:30:41 2019

@author: fergal
"""

from __future__ import print_function
from __future__ import division

import numpy as np
import loadstraws
import pytest

class MonkeyLoadCube(loadstraws.LoadTessCube):
    """A class for testing

    To use the LoadTessCube class normally, we would need access
    to straws from an entire camera/ccd. This is a lot of data
    to haul around entirely for testing purposes.

    This class lets us test all the functionality of LoadTessCube(),
    except for the getStraw() method, which instead of loading a
    file from disk, returns a numpy array of the correct shape.
    This lets me test all the other functionality without the
    need for actual data

    Yes, this is weird and strange. Google Monkey Patching for
    more information.
    """

    def getStraw(self, camera, ccd, col, row):
        shape = (self.nCadences, self.strawSize, self.strawSize)
        straw = np.zeros(shape)
        return straw



def test1():
    """Fixing a bug reported by Susan in testing"""
    obj = MonkeyLoadCube('./testdata')

    cube, col, row = obj.get(1, 1, 221, 250, min_size_pix=40)

    nCad = obj.nCadences
    expect = (nCad, 100, 50)
    assert cube.shape == expect, [cube.shape, expect]



def test_pickABox():

    obj = loadstraws.LoadTessCube(None)
    obj.nCols = 600
    obj.nRows = 700
    obj.nCadences = 40
    obj.strawSize = 50

    bounds = obj.pickBbox(12, 120, 30)
    assert np.all(bounds == np.array([0, 50, 100, 150])), bounds

    bounds = obj.pickBbox(200, 120, 30)
    assert np.all(bounds == np.array([150, 250, 100, 150])), bounds

    bounds = obj.pickBbox(199, 120, 30)
    assert np.all(bounds == np.array([150, 250, 100, 150])), bounds

    bounds = obj.pickBbox(590, 120, 30)
    assert np.all(bounds == np.array([550, 600, 100, 150])), bounds


    bounds = obj.pickBbox(12, 120, 120)
    assert np.all(bounds == np.array([0, 100, 50, 200])), bounds

    bounds = obj.pickBbox(590, 120, 120)
    assert np.all(bounds == np.array([500, 600, 50, 200])), bounds

    with pytest.raises(ValueError):
        bounds = obj.pickBbox(-5, 120, 120)

    with pytest.raises(ValueError):
        bounds = obj.pickBbox(600, 120, 120)

    with pytest.raises(ValueError):
        bounds = obj.pickBbox(60, -120, 120)

    with pytest.raises(ValueError):
        bounds = obj.pickBbox(60, 700, 120)
