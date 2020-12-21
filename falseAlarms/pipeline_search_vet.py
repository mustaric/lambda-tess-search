"""Pipeline code to search and vet TESS data.

"""
import planetSearch as ps
import gen_lightcurve as genlc
import matplotlib.pyplot as plt
import exovetter.tce as TCE
import astropy.units as u
import exovetter.const as const
import lightkurve as lk
from exovetter import vetters


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
                                               tce['duration'].value*24.0,
                                               tce['snr'],
                                               disposition, reason)             
    return st


def vet_all_tces(lc, tce_list, ticid, vetter_list, thresholds, plot=False):
    lcformat = lc.time.format
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
        
        metrics = vet_tce(tce, lc, vetter_list, plot=plot)
        metrics['snr'] = tce['snr']
        disposition, reason = get_disposition(metrics, thresholds)
        result_string = make_result_string(tce, disposition, reason)
        disp_list.append(disposition)
        reason_list.append(reason)
        result_list.append(result_string)
        metrics_list.append(metrics)
        pn=pn+1
        
        
    return result_list, disp_list, reason_list, metrics_list


def search_and_vet_one(ticid, sector, config, vetter_list, plot=True):
    """
    Search and vet one ticid usin gyour config and vetter list
    """
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
    
def plot_lc_tce(tce_list, time, flux, good_time, good_flux, stats):
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

def open_output_file(filename, headerlist, thresholds):
    fobj = open(filename, 'a')
    
    fobj.write("# thresholds: " + str(thresholds)+"\n")
    
    header = headerlist[0]
    for h in headerlist[1:]:
        header = header + ", " + h
        
    fobj.write(header + '\n')
    
    return fobj


