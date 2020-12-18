# -*- coding: utf-8 -*-
"""
Create an eleanor light curve and loop over a BLS
to create TCEs.
Write out the TCEs to a file.
"""

import sys
# Fix path to use local code
sys.path.append('/Users/smullally/Python_Code/lambda-tess-search/planetSearch/')
sys.path.append('/Users/smullally/Python_Code/lambda-tess-search/falseAlarms/')
sys.path[2] = '/Users/smullally/Python_Code/lightkurve/lightkurve'
#%%
import numpy
import planetSearch as ps
import gen_lightcurve as genlc
import matplotlib.pyplot as plt
import exovetter.tce as TCE
import astropy.units as u
import exovetter.const as const
import lightkurve as lk
from exovetter import vetters
from astropy.io import ascii

#%%
def vet_tce(tce, tce_lc, vetter_list, plot=False):

    metrics = dict()
    for v in vetter_list:
        vetter = v
        _ = vetter.run(tce, tce_lc)
        if plot:
            vetter.plot()
        metrics.update(vetter.__dict__)
    return metrics

def get_disposition(metrics, thresholds):
    """Apply thresholds to get a passfail"""
    
    disp = 'PASS' 
    reason = ''
    if metrics['snr'] < thresholds['snr']:
        disp = 'FAIL'
        reason = reason + "-LowSNR-"
    if metrics['norm_lpp'] > thresholds['norm_lpp']:
        disp = 'FAIL'
        reason = reason + "-NormLPP-"
    if metrics['tp_cover'] < thresholds['tp_cover']:
        disp = 'FAIL'
        reason = reason + "-PoorTransitCoverage-"
    if metrics['oe_sigma'] > thresholds['oe_sigma']:
        disp = 'FAIL'
        reason = reason + "-OddEvenDetected-"
    if metrics['sweet']['amp'][0, -1] > thresholds['sweet']:
        disp = 'FAIL'
        reason = reason + "-SWEETHalfPeriod"
    if metrics['sweet']['amp'][1, -1] > thresholds['sweet']:
        disp = 'FAIL'
        reason = reason + "-SWEETAtPeriod"
    if metrics['sweet']['amp'][2,-1] > thresholds['sweet']:
        disp = 'FAIL'
        reason = reason + "-SWEETTwicePeriod-"
    
    
    return disp,reason
    
def make_result_string(tce, disposition, reason):

    st = "%s, %s, %i, %8.4f, %9.4f, %8.3f, %5.3f, %5.2f, %s, %s\n" % \
                                        (tce['target'], tce['event'],
                                               tce['sector'],
                                               tce['period'].value, 
                                               tce['epoch'].value,
                                               tce['depth'].value*1e6,
                                               tce['duration'].value/24.0,
                                               tce['snr'],
                                               disposition, reason)             
    return st


def vet_all_tces(lc, tce_list, vetter_list, plot=False):
    disp_list = []
    reason_list = []
    result_list = []
    metrics_list = []
    pn = 1
    for item in tce_list:
        tce = TCE.Tce(period = item[0]*u.day, epoch=item[1]*u.day, 
                      depth=item[2] * const.frac_amp,
                      duration=item[3]*u.day, 
                      epoch_offset=const.string_to_offset[lcformat],
                      snr=item[4],
                      target = f"TIC {ticid}",
                      sector = lc.sector,
                      event = f"-{pn}-")
        
        metrics = vet_tce(tce, tce_lc, vetter_list, plot=plot)
        metrics['snr'] = tce['snr']
        disposition, reason = get_disposition(metrics, thresholds)
        result_string = make_result_string(tce, disposition, reason)
        disp_list.append(disposition)
        reason_list.append(reason)
        result_list.append(result_string)
        metrics_list.append(metrics)
        
    return result_list, disp_list, reason_list, metrics_list


def search_and_vet_one(ticid, sector, config, vetter_list, plot=True):
    ticid = target['TIC']
    sector = int(target['Sector'])
    lcdata = genlc.hlsp(ticid, sector, author='qlp')
    
    time = lcdata.time.value
    flux = lcdata.flux.value
    flags = lcdata.quality.value

    good_time, meddet_flux = ps.clean_timeseries(time, flux, flags,
                                          config["det_window"], 
                                          config["noise_window"], 
                                          config["n_sigma"], 
                                          sector)
        
    
    tce_list, stats = ps.identifyTces(good_time, meddet_flux, 
                                      bls_durs_hrs=config['bls_dur_hrs'],
                                      minSnr=1, 
                                      fracRemain=0.5, 
                                      maxTces=30, 
                                      minP=config["min_period_days"], 
                                      maxP=config["max_period_days"])
    
    if plot:
        plot_lc_tce(tce_list, time, flux, good_time, meddet_flux)
    
    lcformat = lcdata.time.format
    tce_lc = lk.LightCurve(time=good_time, flux=meddet_flux+1,
                        time_format=lcformat, meta={'sector':sector})
    
    result_strings, disp, reason, metrics_list = vet_all_tces(tce_lc, 
                                                    tce_list, vetter_list,
                                                    plot=plot)
    
    return result_strings, metrics_list
    
def plot_lc_tce(tce_list, time, flux, good_time, good_flux):
    col = ['tab:orange','tab:green','tab:purple','tab:brown',
               'gold','magenta','lightpink']
    plt.figure(figsize=(10,6))
    plt.subplot(211)
    plt.plot(good_time, good_flux,'.', label="clean")
    axes = plt.gca()
    y_min, y_max = axes.get_ylim()
    for n,s in enumerate(stats):
        plt.vlines(stats[n]['transit_times'], y_min, y_max, 
                   colors=col[n], zorder=1)
    
    plt.subplot(212)
    plt.plot(time, flux,'.', label="original")
    
#%%    
#Get List of targets    
target_dir = "/Users/smullally/Science/tess_false_alarms/keplerTargets/target_selection/"
target_file = target_dir + "tic_uncrowded_bright-sectorsObs.csv"

output_dir = "/Users/smullally/Science/tess_false_alarms/vet_results/december2020/"
output_file = output_dir + "disposition.txt"


#%%   
#Set up the cleaningdata prameters and vetter list and vetting thresholds. 

config = {
"det_window" : 45,
"noise_window" : 25,
"n_sigma" : 4,  #noise reject sigma
"max_period_days" : 10,
"min_period_days" : 0.8,
"bls_dur_hrs" : [1,2,4,8,12]
}

vetter_list = [vetters.Lpp(),
                   vetters.OddEven(),
                   vetters.TransitPhaseCoverage(),
                   vetters.Sweet()]

thresholds = {'snr' : 2,
              'norm_lpp' : 2.2,
              'tp_cover' : 0.6,
              'oe_sigma' : 3,
              'sweet' : 3}

#%%
#Here is the loop over different targets
#
target_table = ascii.read(target_file)
fobj = open(output_file, 'a')

choice = (target_table['Sector']==14) | (target_table['Sector']==15)

for target in target_table[choice][0:500]:
    
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
        
        tce_list, stats = ps.identifyTces(good_time, meddet_flux, 
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
    
        result_strings, disp, reason, metrics = vet_all_tces(tce_lc, 
                                                             tce_list, 
                                                             vetter_list)
        
        for r in result_strings:
            print(r)
            fobj.write(r)
    except:
        e = sys.exc_info()[0]
        print(f"Exception raised for TIC {ticid}")
        print(e)
        pass
    
fobj.close()


#%%

ticid = 383724040	#(near to kepler10 but not kepler 10, in S15)
sector = 14
#time, flux, flags = lc.eleanor_corr(ticid, sector)
#Plot a sanity check.
#plt.figure()
#plt.plot(time,flux,'.')

lcdata = genlc.hlsp(ticid, sector, author='qlp')
lcdata.plot()



#%%
#Now work on the BLS

det_window = 29
noise_window = 21
n_sigma = 4  #noise reject sigma

time = lcdata.time.value
flux = lcdata.flux.value
flags = lcdata.quality.value

good_time, meddet_flux = ps.clean_timeseries(time, flux, flags, det_window, \
                                          noise_window, n_sigma, sector)

plt.subplot(211)
plt.plot(good_time, meddet_flux,'.', label="clean")
plt.subplot(212)
plt.plot(time,flux,'.', label="original")

#%%    
tce_list, stats = ps.identifyTces(good_time, meddet_flux, bls_durs_hrs=[1,2,4,8,12],
                                 minSnr=1, fracRemain=0.5, 
                                 maxTces=30, minP=None, maxP=None)
#%%

vetter_list = [vetters.Lpp(),
                   vetters.OddEven(),
                   vetters.TransitPhaseCoverage(),
                   vetters.Sweet()]

thresholds = {'snr' : 3,
              'norm_lpp' : 1,
              'tp_cover' : 0.5,
              'oe_sigma' : 3,
              'sweet' : 2.5}
    


lcformat = lcdata.time.format
tce_lc = lk.LightCurve(time=good_time, flux=meddet_flux+1,
                        time_format=lcformat)

result_strings, disp, reason = vet_all_tces(tce_lc, tce_list, vetter_list)    

#%%
#Write the output for this one light curve.
output_dir = "/Users/smullally/Science/tess_false_alarms/vet_results/december2020/"
output_file = output_dir + "disposition.txt"
fobj = open(output_file, 'a')
for r in result_strings:
    print(r)
    fobj.write(r)
    
fobj.close()

#%%
#Output
import pickle
import numpy as np

ddir = "/Users/smullally/TESS/planetSearch/keplerFA/test1/"
resultfile = "%s/tic%014u-s%02u-result.txt" %(ddir,ticid,sector)
statsfile = "%s/tic%014u-s%02u-stats.pkl" % (ddir,ticid, sector)
print(resultfile, statsfile)

with open(statsfile, 'w') as outfile:
    pickle.dump(stats, outfile)

np.savetxt(resultfile, results)

