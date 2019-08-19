import json
import numpy as np
import requests
import warnings
from astropy.utils.data import CacheMissingWarning
from astropy.config import ConfigurationMissingWarning
from astropy.config import InvalidConfigurationItemWarning
warnings.simplefilter('ignore', CacheMissingWarning)
warnings.simplefilter('ignore', ConfigurationMissingWarning)
warnings.simplefilter('ignore', InvalidConfigurationItemWarning)
from astropy.io import fits

import boto3
import os
#from dotenv import load_dotenv


#from astroquery.mast import Observations

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
    afile = event['fileloc']
    head = fits.getheader(afile, 1)
    
    midtime = head['TSTART'] + (head['TSTOP'] - head['TSTART'])/2
    
    print(midtime)
    
    filename="/tmp/myfile.fits"
    #url="tess/public/mast/tess-s0001-1-1-cube.fits"
    #bucket = s3.Bucket('stpubdata')
    url="tess2018206192942-s0001-1-1-0120-s_ffic.fits"
    s3 = boto3.resource('s3')

    bucket = s3.Bucket('firstlambdacaom')
    
    lightcurve = np.zeros(10)
    for i in range(10):
        bucket.download_file(url, filename, ExtraArgs={"RequestPayer": "requester"})
        data = fits.getdata(filename)
        point = np.nanmedian(data[100:150,100:150])
        lightcurve[i] = point
        os.remove(filename)
    
    out_file="/tmp/output.csv"
    np.savetxt(out_file,lightcurve, delimiter=',')
    s3_client = boto3.client('s3')
    resp=s3_client.upload_file(out_file, "firstlambdacaom", "output.csv")
    
    return {
            "statusCode": 200,
            "body": json.dumps({
                    "midtime": midtime,
                    "response": str(resp),
                    "exposure": head['EXPOSURE'],
                    "lightcurve": str(list(lightcurve))
                    
                    })
            }
    
    
    
    """    
    obsTable = Observations.query_object(name, radius = ".005 deg")
    
    if len(obsTable) > 0:
        filter2 = obsTable['target_name'] == 'TESS FFI'
        print("Tess-Sector-Camera-CCD for %s" % name)
        print(obsTable[filter2]['obs_id'])
        table = obsTable[filter2]['obs_id']
    else:
        table = "None"

    return {
        "statusCode": 200,
        "body": json.dumps({
            "message": "hello TIC %s !" % tic 
          #  "ffis": ", ".join(table)
        })
    }
    
    """    
#event = {'fileloc': "https://archive.stsci.edu/missions/tess/ffi/s0002/2018/235/1-2/tess2018235152941-s0002-1-2-0121-s_ffic.fits"}
#context ={}
#retvalue=lambda_handler(event, context)
