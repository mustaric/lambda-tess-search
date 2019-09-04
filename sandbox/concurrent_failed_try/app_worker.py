import json
import os

import boto3
import numpy as np
from astropy.io import fits
from astropy.utils.misc import JsonCustomEncoder


# https://stackoverflow.com/questions/49330080/numpy-2d-array-selecting-indices-in-a-circle
def circular_mask(arr_shape, cx, cy, r):
    y = np.arange(arr_shape[0])
    x = np.arange(arr_shape[1])
    return ((x[np.newaxis, :] - cx) ** 2 + (y[:, np.newaxis] - cy) ** 2) < (r ** 2)  # noqa


def lambda_handler(event, context):
    """Extract light curve data from one TESS full frame image.

    Parameters
    ----------
    event : dict
        API Gateway Lambda Proxy Input Format.
        Event doc: https://docs.aws.amazon.com/apigateway/latest/developerguide/set-up-lambda-proxy-integrations.html#api-gateway-simple-proxy-for-lambda-input-format

    context : object
        Lambda Context runtime methods and attributes.
        Context doc: https://docs.aws.amazon.com/lambda/latest/dg/python-context-object.html

    Returns
    ------
    result : dict
        API Gateway Lambda Proxy Output Format.
        Return doc: https://docs.aws.amazon.com/apigateway/latest/developerguide/set-up-lambda-proxy-integrations.html

    """  # noqa
    # Example:
    # tess/public/ffi/s0001/2018/206/1-1/tess2018206192942-s0001-1-1-0120-s_ffic.fits
    key = event['key']
    # For signal
    xpos = int(event['xpos'])
    ypos = int(event['ypos'])
    radius = int(event['radius'])
    # For background
    threshold = float(event['bright_pixel_threshold'])

    filename = os.path.join('/tmp', key.split('/')[-1])

    s3 = boto3.resource('s3')
    bucket = s3.Bucket('stpubdata')
    bucket.download_file(
        key, filename, ExtraArgs={"RequestPayer": "requester"})

    np.seterr(all='ignore')

    with fits.open(filename) as pf:
        # EXT 0: Primary header
        obstime = 0.5 * (pf[0].header['TSTOP'] + pf[0].header['TSTART'])  # BJD

        # EXT 1: cal
        # Signal is inside a circle.
        mask = circular_mask(pf[1].data.shape, xpos, ypos, radius)
        signal = np.nanmedian(pf[1].data[mask])
        # Background is everything else.
        bg_arr = pf[1].data[~mask]
        background = np.nanmedian(bg_arr[(bg_arr * threshold) < signal])

    # Clean up
    os.remove(filename)

    # signal and background would be NaN if no valid data found, unless
    # exception is raised, then Lambda would return something else, not this.
    return {"statusCode": 200,
            "body": json.dumps({
                "midtime": obstime,
                "signal": signal,
                "background": background}, cls=JsonCustomEncoder)}
