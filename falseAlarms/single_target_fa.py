# -*- coding: utf-8 -*-
"""
Create an eleanor light curve and loop over a BLS
to create TCEs.
Write out the TCEs to a file.
"""

import numpy
import planetSearch as ps
import gen_lightcurve as lc
import matplotlib.pyplot as plt

ticid = 350814695	#(near to kepler10 but not kepler 10, in S15)
sector = 15
time, flux, flags = lc.eleanor_corr(ticid, sector)

#%%
#Plot a sanity check.
plt.figure()
plt.plot(time,flux,'.')



#%%
#Now work on the BLS

det_window = 29
noise_window = 21
n_sigma = 4  #noise reject sigma

good_time, meddet_flux = ps.clean_timeseries(time, flux, flags, det_window, \
                                          noise_window, n_sigma, sector)

plt.subplot(211)
plt.plot()
plt.plot(good_time, meddet_flux,'.')
plt.subplot(212)
plt.plot(time,flux,'.')

#%%    
    
results, stats = ps.identifyTces(good_time, meddet_flux, bls_durs_hrs=[1,2,4,8,12],
                                 minSnr=1, fracRemain=0.5, 
                                 maxTces=30, minP=None, maxP=None)
#result contains bls_period, bls_t0, bls_depth, bls_duration, simple_snr
print(results)

#%%
#Output
import pickle

ddir = "/Users/smullally/TESS/planetSearch/keplerFA/test1/"
resultfile = "%s/tic%014u-s%02u-result.txt" %(ddir,ticid,sector)
statsfile = "%s/tic%014u-s%02u-stats.pkl" % (ddir,ticid, sector)
print(resultfile, statsfile)

with open(statsfile, 'w') as outfile:
    pickle.dump(stats, outfile)

np.savetxt(resultfile,results)

