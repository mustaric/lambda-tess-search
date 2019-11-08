import json

import boto3

import planetSearch as ps
import data_io as io
import os


def lambda_handler(event, context):
    """Sample pure Lambda function

    Parameters
    ----------
    event: dict, required
        API Gateway Lambda Proxy Input Format

        Event doc: https://docs.aws.amazon.com/apigateway/latest/developerguide/set-up-lambda-proxy-integrations.html#api-gateway-simple-proxy-for-lambda-input-format

    context: object, required
        Lambda Context runtime methods and attributes

        Context doc: https://docs.aws.amazon.com/lambda/latest/dg/python-context-object.html

    Returns
    ------
    API Gateway Lambda Proxy Output Format: dict

        Return doc: https://docs.aws.amazon.com/apigateway/latest/developerguide/set-up-lambda-proxy-integrations.html
    """
    print(event)
  
    #Input variables
    aminP = 0.75
    amaxP = 10.0
    durs = [1, 2, 3, 4, 5, 6]
    min_snr = 4
    max_tce = 3
    frac_remain = 0.8
    det_window = 45
    noise_window = 12
    n_sigma = 4
    search_bucket = "tesssearchresults"
    #---------
    
    cloud = True
    #Local Storage
    local_filename = "/tmp/mylightcurve.fits"
    local_detrend_fn = "/tmp/detrended.fits"
    out_file = "/tmp/output.csv"
    
    b_filename = event['Records'][0]['s3']['object']['key']
    bucket = event['Records'][0]['s3']['bucket']['name']   
    
    #If not in the cloud bucket begins with /, do the following to set up for a test.
    if bucket[0] == "/":
        cloud = False  #For testing
        local_filename = bucket + b_filename
        local_dir = "/Users/smullally/TESS/lambdaSearch/test/tesssearchresults/"
        local_detrend_fn = local_dir + "detrended.fits"
        out_file = local_dir + "outfile.fits"
    
    meta = dict()
    
    #Get the information from the light curve file.
    if cloud:
        time, flux, qflags, phead = io.read_lightcurve_lambda(bucket, \
                                                    b_filename, local_filename)
    else:
        time, flux, qflags, phead = io.read_lightcurve_lambda_local(bucket, \
                                            b_filename, local_filename)
    
    print(time,flux,qflags)
    ticid, camera, sector, ccd = io.read_header(phead)

    namestr = "tic%012u/tic%012u_s%04u-%1u-%1u" % \
            (int(ticid), int(ticid), int(sector),int(camera), int(ccd))
    
    #Detrend
    good_time, meddet_flux = ps.clean_timeseries(time, flux, qflags, det_window, \
                                          noise_window, n_sigma)
    
    #Take BLS
    results, stats = ps.identifyTces(good_time, meddet_flux, bls_durs_hrs=durs,\
                                     minSnr=min_snr, \
                                     fracRemain=frac_remain,\
                                     maxTces=max_tce, minP=aminP, maxP=amaxP)
    print(results)
    
    #Now write out results.
    
    bucket_out_name = namestr + "_plsearch" + '.csv'
    bucket_detrend_name = namestr + "_detrend" + '.fits'
    
    
    io.write_results(out_file, int(ticid), results, stats, **meta)
    io.write_timeseries(local_detrend_fn, good_time, meddet_flux, phead)
    
    if cloud:  
        #Write to the S3 bucket.
        s3_client = boto3.client('s3')
        resp = s3_client.upload_file(out_file, search_bucket, bucket_out_name)
        resp = s3_client.upload_file(local_detrend_fn, search_bucket, bucket_detrend_name)
    else:
        resp = "not cloud"
        if not os.path.exists(local_dir + "tic%012u" % (int(ticid))):
            os.mkdir(local_dir + "tic%012u" % (int(ticid)))
        try:
            os.remove(local_dir + bucket_out_name)
        except:
            pass
        os.rename(out_file, local_dir + bucket_out_name)
        try:
            os.remove(local_dir + bucket_detrend_name)
        except:
            pass
        os.rename(local_detrend_fn, local_dir + bucket_detrend_name)
    
    return {
            "statusCode": 200,
            "body": json.dumps({
                    "outname": bucket_out_name,
                    "response": str(resp),
                    "period": str(results[0][0]),
                    "epoch": str(results[0][1])
                    })
            }
  
    
def test1():
    
    event = {
    "Records": [
    {
      "s3": {
        "s3SchemaVersion": "1.0",
        "configurationId": "b0efd5b1-cc92-47b4-8501-1c34f5eba235",
        "bucket": {
          "name": "/Users/smullally/TESS/lambdaSearch/test/ffi-lightcurves/"
        },
        "object": {
          "key": "tic000025155310_s0001-4-1_lcc.fits"
        }
      }
    }
  ]
}
    context = {}
    
    out = lambda_handler(event, context)
    
    assert out["statusCode"] == 200


def test2():
    
    event = {
    "Records": [
    {
      "s3": {
        "s3SchemaVersion": "1.0",
        "configurationId": "b0efd5b1-cc92-47b4-8501-1c34f5eba235",
        "bucket": {
          "name": "/tmp/"
        },
        "object": {
          "key": "tic000147424426_s0001-1-1_stlc.fits"
        }
      }
    }
  ]
}
    context = {}
    
    out = lambda_handler(event, context)
    
    assert out["statusCode"] == 200