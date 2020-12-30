# -*- coding: utf-8 -*-
"""
Generate TCEs from HLSP light curves 
Vet those TCEs
Do this for a list of TIC IDs that overlap Kepler Targets.
Write out to file.
"""

import sys
# Fix path to use local code
sys.path.append('/Users/smullally/Python_Code/lambda-tess-search/planetSearch/')
sys.path.append('/Users/smullally/Python_Code/lambda-tess-search/falseAlarms/')
sys.path[2] = '/Users/smullally/Python_Code/lightkurve/lightkurve'

import planetSearch as ps
import gen_lightcurve as genlc
import matplotlib.pyplot as plt
import lightkurve as lk
from exovetter import vetters
import exovetter.tce as TCE
from astropy.io import ascii
import pipeline_search_vet as pipe
   
#%%    
#Get List of targets    
target_dir = "/Users/smullally/Science/tess_false_alarms/keplerTargets/target_selection/"
target_file = target_dir + "tic_bright14-sectorsObs-all.csv"

output_dir = "/Users/smullally/Science/tess_false_alarms/vet_results/search_results_2020dec30/"
output_file = output_dir + "disposition_tic_bright14-sectorsObs-2.csv"

run_descriptor = "2020dec30"

#List of output information
headerlist = ['target', 'pn', 'sector', 'period', 'epoch', 'depth', 
              'duration', 'snr', 'dispotion', 'reason']

#%%   
#Set up the cleaningdata prameters and vetter list and vetting thresholds. 

config = {
"det_window" : 65,
"noise_window" : 27,
"n_sigma" : 4.5,  #noise reject sigma
"max_period_days" : 10,
"min_period_days" : 0.8,
"bls_dur_hrs" : [1,2,4,8,12]
}

vetter_list = [vetters.Lpp(),
                   vetters.OddEven(),
                   vetters.TransitPhaseCoverage(),
                   vetters.Sweet()]

thresholds = {'snr' : 1,
              'norm_lpp' : 2.0,
              'tp_cover' : 0.6,
              'oe_sigma' : 3,
              'sweet' : 3}
 
target_table = ascii.read(target_file)
fobj = pipe.open_output_file(output_file, headerlist, thresholds)
fobj.close()

#%%
#Now loop over targets in your list.
#choice = (target_table['Sector']==14) | (target_table['Sector']==15)
choice = (target_table['Sector']==14)
for target in target_table[choice][0:200]:
    
    ticid = target['TIC']
    sector = int(target['Sector'])
    try:
        lcdata = genlc.hlsp(ticid, sector, author='qlp')
        
        time = lcdata.time.value
        flux = lcdata.flux.value
        flags = lcdata.quality.value
    
        good_time, meddet_flux = ps.clean_timeseries(time, flux, flags,
                                              config["det_window"], 
                                              config["noise_window"], 
                                              config["n_sigma"], 
                                              sector)
        
        tce_dlist, stats = ps.identifyTces(good_time, meddet_flux, 
                                          bls_durs_hrs=config['bls_dur_hrs'],
                                          minSnr=1, 
                                          fracRemain=0.5, 
                                          maxTces=30, 
                                          minP=config["min_period_days"], 
                                          maxP=config["max_period_days"])
        
        lcformat = lcdata.time.format
        tce_lc = lk.LightCurve(time=good_time, flux=meddet_flux+1,
                            time_format=lcformat, meta={'sector':sector})
        #tce_lc['sector'] = sector
    
        result_strings, disp, reason, metrics, tces = pipe.vet_all_tces(tce_lc, 
                                                             tce_dlist, ticid,
                                                             vetter_list,
                                                             thresholds)
        for r in result_strings:
            print(r)
            fobj = open(output_file, 'a')
            fobj.write(r)
            fobj.close()
        for tce in tces:
            tcefilename = "tic%09i-%02i-%s.json" % (int(tce['target'][5:]), 
                                                    int(tce['event']), 
                                                    run_descriptor)
            print(tcefilename)
            full_filename = output_dir + tcefilename
            tce.to_json(full_filename)
    except:
        e = sys.exc_info()[0]
        print(f"Exception raised for TIC {ticid}")
        print(e)
        pass
    
fobj.close()


#%%

# ticid = 383724040	#(near to kepler10 but not kepler 10, in S15)
# sector = 14
# #time, flux, flags = lc.eleanor_corr(ticid, sector)
# #Plot a sanity check.
# #plt.figure()
# #plt.plot(time,flux,'.')

# lcdata = genlc.hlsp(ticid, sector, author='qlp')
# lcdata.plot()



#%%
#Now work on the BLS

# det_window = 29
# noise_window = 21
# n_sigma = 4  #noise reject sigma

# time = lcdata.time.value
# flux = lcdata.flux.value
# flags = lcdata.quality.value

# good_time, meddet_flux = ps.clean_timeseries(time, flux, flags, det_window, \
#                                           noise_window, n_sigma, sector)

# plt.subplot(211)
# plt.plot(good_time, meddet_flux,'.', label="clean")
# plt.subplot(212)
# plt.plot(time,flux,'.', label="original")

# #%%    
# tce_list, stats = ps.identifyTces(good_time, meddet_flux, bls_durs_hrs=[1,2,4,8,12],
#                                  minSnr=1, fracRemain=0.5, 
#                                  maxTces=30, minP=None, maxP=None)
# #%%

# vetter_list = [vetters.Lpp(),
#                    vetters.OddEven(),
#                    vetters.TransitPhaseCoverage(),
#                    vetters.Sweet()]

# thresholds = {'snr' : 3,
#               'norm_lpp' : 1,
#               'tp_cover' : 0.5,
#               'oe_sigma' : 3,
#               'sweet' : 2.5}
    


# lcformat = lcdata.time.format
# tce_lc = lk.LightCurve(time=good_time, flux=meddet_flux+1,
#                         time_format=lcformat)

# result_strings, disp, reason = vet_all_tces(tce_lc, tce_list, vetter_list)    

# #%%
# #Write the output for this one light curve.
# output_dir = "/Users/smullally/Science/tess_false_alarms/vet_results/december2020/"
# output_file = output_dir + "disposition.txt"
# fobj = open(output_file, 'a')
# for r in result_strings:
#     print(r)
#     fobj.write(r)
    
# fobj.close()

# #%%
# #Output
# import pickle
# import numpy as np

# ddir = "/Users/smullally/TESS/planetSearch/keplerFA/test1/"
# resultfile = "%s/tic%014u-s%02u-result.txt" %(ddir,ticid,sector)
# statsfile = "%s/tic%014u-s%02u-stats.pkl" % (ddir,ticid, sector)
# print(resultfile, statsfile)

# with open(statsfile, 'w') as outfile:
#     pickle.dump(stats, outfile)

# np.savetxt(resultfile, results)

