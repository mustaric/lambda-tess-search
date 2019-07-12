"""Memory profiling for partial cube retrieval from S3 on EC2.

This is basically profiling the example that is given in the
docstring of ``aws_remote_fits.py``.

To run::

    python measure_memory_aws.py

Line #    Mem usage    Increment   Line Contents
================================================
    14     59.1 MiB     59.1 MiB   @profile
    15                             def s3_subcube():
    16     59.1 MiB      0.0 MiB       s3_client = boto3.client(
    17     59.1 MiB      0.0 MiB           's3',
    18     59.1 MiB      0.0 MiB           aws_access_key_id=...,
    19     63.5 MiB      4.4 MiB           aws_secret_access_key=...)
    20     63.5 MiB      0.0 MiB       s3_file = RemoteAWSFITSFile(
    21     64.9 MiB      1.4 MiB           s3_client, 'stpubdata', 'tess/public/mast/tess-s0001-1-1-cube.fits')
    22     65.4 MiB      0.5 MiB       subcube = s3_file[slice(0, 10), slice(0,10), slice(0, 1), slice(0, 1)]
    23     65.4 MiB      0.0 MiB       subcube = subcube.reshape(10, 10).T
    24     65.6 MiB      0.2 MiB       x = subcube.mean()  # noqa
    25     65.7 MiB      0.0 MiB       fits.writeto('subcube.fits', subcube, overwrite=True)
    26     65.7 MiB      0.0 MiB       with fits.open('subcube.fits') as pf:
    27     65.7 MiB      0.0 MiB           data = pf[0].data.copy()
    28     65.7 MiB      0.0 MiB       np.testing.assert_array_equal(subcube, data)

Line #    Mem usage    Increment   Line Contents
================================================
    14     59.1 MiB     59.1 MiB   @profile
    15                             def s3_subcube():
    16     59.1 MiB      0.0 MiB       s3_client = boto3.client(
    17     59.1 MiB      0.0 MiB           's3',
    18     59.1 MiB      0.0 MiB           aws_access_key_id=...,
    19     63.3 MiB      4.3 MiB           aws_secret_access_key=...)
    20     63.3 MiB      0.0 MiB       s3_file = RemoteAWSFITSFile(
    21     64.7 MiB      1.4 MiB           s3_client, 'stpubdata', 'tess/public/mast/tess-s0001-1-1-cube.fits')
    22     66.0 MiB      1.3 MiB       subcube = s3_file[slice(0, 100), slice(0, 100), slice(0, 1), slice(0, 1)]
    23     66.0 MiB      0.0 MiB       subcube = subcube.reshape(100, 100).T
    24     66.0 MiB      0.0 MiB       x = subcube.mean()  # noqa
    25     66.1 MiB      0.1 MiB       fits.writeto('subcube.fits', subcube, overwrite=True)
    26     66.1 MiB      0.0 MiB       with fits.open('subcube.fits') as pf:
    27     66.1 MiB      0.0 MiB           data = pf[0].data.copy()
    28     66.1 MiB      0.0 MiB       np.testing.assert_array_equal(subcube, data)

For completeness, here is the memory footprint after running
the two memory profiling above on EC2 with free-tier RHEL::

    $ free -h
                  total        used        free      shared  buff/cache   available
    Mem:          819Mi       164Mi       110Mi        10Mi       544Mi       508Mi
    Swap:            0B          0B          0B

In the end, I think memory is of no concern as run time becomes
a bottleneck when sub-cube becomes larger than 10x10.

"""
# STDLIB
import os

# THIRD-PARTY
import boto3
import numpy as np
from astropy.io import fits
from memory_profiler import profile

# LOCAL
from aws_remote_fits import RemoteAWSFITSFile


@profile
def s3_subcube():
    s3_client = boto3.client(
        's3',
        aws_access_key_id=os.environ['aws_access_key_id'],
        aws_secret_access_key=os.environ['aws_secret_access_key'])
    s3_file = RemoteAWSFITSFile(
        s3_client, 'stpubdata', 'tess/public/mast/tess-s0001-1-1-cube.fits')
    subcube = s3_file[slice(0, 10), slice(0, 10), slice(0, 1), slice(0, 1)]
    subcube = subcube.reshape(10, 10).T
    x = subcube.mean()  # noqa
    fits.writeto('subcube.fits', subcube, overwrite=True)
    with fits.open('subcube.fits') as pf:
        data = pf[0].data.copy()
    np.testing.assert_array_equal(subcube, data)


if __name__ == '__main__':
    s3_subcube()
