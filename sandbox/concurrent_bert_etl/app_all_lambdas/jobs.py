import boto3
import os
import typing

import numpy as np

from astropy.io import fits
from astroquery.mast import Observations

from bert import binding, constants, utils

#from tess_bert import shortcuts

from urllib.parse import urlparse, ParseResult

OBS_ID: str = 'tess-s0001-1-1'
DATA_DIR: str = os.path.join(os.getcwd(), 'data')
if not os.path.exists(DATA_DIR):
    os.makedirs(DATA_DIR)

@binding.follow('noop')
def aws_fullframe_fits():
    """
    Loop through full frame files, extract a subarray, and calculate mean.
    This must be done in a way that the file is deleted as soon as it is, no longer
    necessary to keep, so we do not use up all the disk space.
    """
    work_queue, done_queue, ologger = utils.comm_binders(aws_fullframe_fits)

    obs_table = Observations.query_criteria(obs_id=OBS_ID)
    products = Observations.get_product_list(obs_table)
    filtered = Observations.filter_products(products, productSubGroupDescription="FFIC", mrp_only=False)
    Observations.enable_cloud_dataset()

    for idx, s3_url in enumerate(Observations.get_cloud_uris(filtered, include_bucket=True)):
        url_parts: ParseResult = urlparse(s3_url)
        filepath: str = os.path.join(DATA_DIR, os.path.basename(url_parts.path))
        done_queue.put({
            'bucket_path': url_parts.path.strip('/'),
            'filepath': filepath,
            'bucket': url_parts.netloc
        })
        if idx > 2 and constants.DEBUG:
            break

@binding.follow(aws_fullframe_fits, pipeline_type=constants.PipelineType.CONCURRENT)
def find_mean():
    work_queue, done_queue, ologger = utils.comm_binders(find_mean)
    for details in work_queue:
        bucket: 'boto3.resources.factory.s3.Bucket' = boto3.resource('s3').Bucket(name=details['bucket'])
        filename: str = os.path.basename(details['filepath'])
        ologger.info(f'Downloading File[{filename}]')
        bucket.download_file(details['bucket_path'], details['filepath'], {'RequestPayer': 'requester'})
        with fits.open(details['filepath']) as pf:
            details['indiv_mean'] = str(pf[1].data[:10, :10].mean())
            done_queue.put(details)
            # details['midtime'] = 0.5 *(pf[0].header['TSTOP'] + pf[0].header['TSTART'])
            # mask = shortcuts.circular_mask(pf[1].data.shape, xpos, ypos, radius)


        os.remove(details['filepath'])

@binding.follow(find_mean)
def find_final_mean():
    work_queue, done_queue, ologger = utils.comm_binders(find_final_mean)
    indiv_means: typing.List[float] = []
    for details in work_queue:
        indiv_means.append(details['indiv_mean'])

    final_mean: typing.Any = np.mean([item for item in map(lambda x: float(x), indiv_means)])
    ologger.info(f'Final Mean[{final_mean}]')
