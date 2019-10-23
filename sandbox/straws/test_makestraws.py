# -*- coding: utf-8 -*-
# Copyright 2017-2018 Orbital Insight Inc., all rights reserved.
# Contains confidential and trade secret information.
# Government Users: Commercial Computer Software - Use governed by
# terms of Orbital Insight commercial license agreement.

"""
Created on Wed Aug 21 10:02:05 2019

@author: fergal
"""

from __future__ import print_function
from __future__ import division

from pdb import set_trace as debug
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import matplotlib as mpl
import pandas as pd
import numpy as np

import makestraws
import pytest

def test_pickABox():

    obj = makestraws.LoadTessStraw(None)
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
