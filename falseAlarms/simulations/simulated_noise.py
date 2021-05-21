#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 12:29:05 2021

@author: smullally

This reads in the simulated TCEs (gaussian noise and searched by corazon bls)
Then it plots histogram of those TCEs snr.

"""
import pandas as p
import numpy as np
import matplotlib.pyplot as plt

sim_datafiles = ["/Users/smullally/Science/tess_false_alarms/simulations/simulated_tces_200_1200_20210215_v1.txt",
                 "/Users/smullally/Science/tess_false_alarms/simulations/simulated_tces_200_1200_20210215_v2.txt",
                 "/Users/smullally/Science/tess_false_alarms/simulations/simulated_tces_200_400_20210219_v3.txt"]

names = ['ID', 'event', 'sector', 'period', 'epoch', 'depth', 'duration', 'snr', 'disp', 'reason', 'noise']

sim = p.DataFrame(columns=names)
for file in sim_datafiles:
    simdf = p.read_csv(file, names=names, skiprows=1)
    sim = sim.append(simdf)
    
#sim1 = p.read_csv(sim_datafile, names=names, skiprows=1)
#sim2 = p.read_csv(sim_datafile, names=names, skiprows=1)

fps = sim['disp'] == 'FAIL'
lownoise = sim['noise'] < np.mean(sim['noise'])
highnoise = ~lownoise

nbins = np.arange(0,1.5, .1)
plt.figure()
plt.subplot(311)
plt.hist(sim['snr'], bins=nbins, label="All")
plt.legend()

plt.subplot(312)
plt.hist(sim[lownoise]['snr'], bins=nbins)
plt.subplot(313)
plt.hist(sim[highnoise]['snr'], bins=nbins, label = ['high noise'])
plt.legend()
