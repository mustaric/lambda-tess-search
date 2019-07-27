"""Adapted from https://mast-labs.stsci.io/2018/12/tess-data-available-on-aws .
This explores the feasibility of extracting data from a series of
full frame images instead of from a large cube.

For this to work, ``~/.aws/credentials`` is needed in the EC2 instance.

To run::

    python aws_fullframe_fits.py

On a RHEL free-tier EC2 instance, this gave a subcube mean of 0.010893357,
which is meaningless except to show that it executed as expected.
The function took a whooping 5146.88 seconds, which is almost 1.5 hours,
to process one observation ID. Disk space used was about 1.5 GB (out of 10 GB
available), which was taken up mostly by installed software. Memory usage
was well within the allocated limit (few hundred MB).

So, what does this mean for AWS Lambda? High-volume S3 download does not seem
practical. Mounting it as FUSE is not possible on Lambda (only EC2) and
will not magically turn it into Central Storage anyway.

"""
import os
import time

import boto3
import numpy as np
from astropy.io import fits
from astroquery.mast import Observations

# NOTE: Use your own key values here.
os.environ['AWS_ACCESS_KEY_ID'] = 'somekey'
os.environ['AWS_SECRET_ACCESS_KEY'] = 'somesecret'

# NOTE: Change TESS observation ID as needed.
obs_id = 'tess-s0001-1-1'

# Find full frame dataset for the observation ID.
obs_table = Observations.query_criteria(obs_id=obs_id)
products = Observations.get_product_list(obs_table)
filtered = Observations.filter_products(
    products, productSubGroupDescription="FFIC", mrp_only=False)

# Set up AWS S3 bucket to pull data from.
Observations.enable_cloud_dataset()
s3_urls = Observations.get_cloud_uris(filtered, include_bucket=False)
s3 = boto3.resource('s3')
bucket = s3.Bucket('stpubdata')


def time_mean():
    """Loop through full frame files, extract a subarray, and calculate mean.
    This must be done in a way that the file is deleted as soon as it is
    no longer necessary to keep, so we do not use up all the disk space.

    .. note:: Algorithm can also be modified to construct subarrays
              into a subcube.

    Returns
    -------
    final_mean : float
        Mean of the extracted subcube.

    total_time : float
        Run time in seconds.

    """
    indiv_mean = []
    extra_args = {"RequestPayer": "requester"}

    t_start = time.time()  # Start timer

    # NOTE: Over 1000 files total, use indexing for smaller loop.
    for url in s3_urls:
        filename = url.split('/')[-1]
        bucket.download_file(url, filename, ExtraArgs=extra_args)

        with fits.open(filename) as pf:
            # EXT 1: cal
            # NOTE: Modify subarray to extract as needed.
            x = pf[1].data[:10, :10].mean()
            indiv_mean.append(x)

        os.remove(filename)

    final_mean = np.mean(indiv_mean)

    # Stop timer and calculate run time.
    t_end = time.time()
    total_time = t_end - t_start  # seconds

    return final_mean, total_time


if __name__ == '__main__':
    final_mean, total_time = time_mean()
    print(f'Data mean: {final_mean}\nRun time: {total_time} s')
