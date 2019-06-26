"""Adapted from code written by Tom Robitaille (@astrofrog) at
https://gist.github.com/astrofrog/404a5682311bae221fd5308c1498da94

Requires:

* boto3
* numpy
* astropy

Examples
--------
Import all the necessary modules:

>>> import boto3
>>> import numpy as np
>>> from astropy.io import fits

This imports the data class from a local module named
``aws_remote_fits.py``. If you have copied the code
for ``RemoteAWSFITSFile`` directly into the active Python
session, you can skip this:

>>> from aws_remote_fits import RemoteAWSFITSFile

Authenticate with AWS S3. Replace ``'xxx'`` with
your credentials:

>>> s3_client = boto3.client(
...     's3', aws_access_key_id='xxx', aws_secret_access_key='xxx')

Create an instance of ``RemoteAWSFITSFile`` that holds
the desired TESS cube file being hosted on S3:

>>> s3_file = RemoteAWSFITSFile(
...     s3_client, 'stpubdata', 'tess/public/mast/tess-s0001-1-1-cube.fits')

On a successful read, you will see some attributes being
populated:

>>> print(s3_file.header)
TENSION= 'IMAGE   '           / Image extension
BITPIX  =                  -32 / array data type
NAXIS   =                    4 / number of array dimensions
NAXIS1  =                    2
NAXIS2  =                 1282
NAXIS3  =                 2136
NAXIS4  =                 2078
PCOUNT  =                    0 / number of parameters
GCOUNT  =                    1 / number of groups
END
...
>>> print(s3_file.shape)
(2078, 2136, 1282, 2)

This is the part where ranged GET request to S3 bucket happens.
This is almost the same as ``cube[:10, :10, 0, 0]``, except that
the result returned needs a reshape and a transpose afterwards:

>>> subcube = s3_file[slice(0, 10), slice(0, 10), slice(0, 1), slice(0, 1)]
>>> subcube = subcube.reshape(10, 10).T

Once you have the sub-cube, you can perform regular Numpy
operations on it:

>>> print(subcube.mean())
0.009653374

To write it to a file locally as a simple single-extension FITS:

>>> fits.writeto('subcube.fits', subcube, overwrite=True)

To make sure it round-trips:

>>> with fits.open('subcube.fits') as pf:
...     data = pf[0].data.copy()
>>> np.testing.assert_array_equal(subcube, data)

"""
import numpy as np
from astropy.io import fits


class RemoteAWSFITSFile:
    """Class to handle partial download of FITS file
    containing data cube stored in AWS S3 bucket.

    Currently, this only handles the following file format, where
    ``ImageHDU`` has 4D array::

        No.    Name      Ver    Type      Cards   Dimensions   Format
        0  PRIMARY       1 PrimaryHDU      31     ()
        1                1 ImageHDU         9     (..., ...)   float32

    Parameters
    ----------
    s3_client : obj
        An instance of authenticated ``boto3.client``.

    s3_bucket : str
        S3 bucket name. Example: ``'stpubdata'``

    s3_key : str
        S3 key to the file.
        Example: ``'tess/public/mast/tess-s0001-1-1-cube.fits'``

    """
    def __init__(self, s3_client, s3_bucket, s3_key):

        self.s3_client = s3_client
        self.s3_bucket = s3_bucket
        self.s3_key = s3_key

        # Load complete header which allows us to find start of data

        header_block_size = 2879

        # EXT=0 with primary header only, no data.
        hdr_str = self.get_data_single_range(0, header_block_size)
        self.primary_header = fits.Header.fromstring(hdr_str.decode('ascii'))

        # EXT=1 with header and data, other extensions ignored.
        hdr_str = self.get_data_single_range(
            header_block_size + 1, 2 * header_block_size + 1)
        self.header = fits.Header.fromstring(hdr_str.decode('ascii'))

        self.data_start = 2 * header_block_size + 2

        # Data type, see
        # https://docs.astropy.org/en/stable/io/fits/usage/image.html
        bitpix = self.header['BITPIX']
        if bitpix == 8:
            data_type = np.uint8
        elif bitpix == 16:
            data_type = np.int16
        elif bitpix == 32:
            data_type = np.int32
        elif bitpix == -32:
            data_type = np.float32
        elif bitpix == -64:
            data_type = np.float64
        else:
            raise ValueError('BITPIX={} not supported'.format(bitpix))

        self.data_type = np.dtype(data_type).newbyteorder('>')
        self.itemsize = np.dtype(data_type).itemsize

        # When fitsinfo gives (2, 1282, 2136, 2078),
        # shape gives (2078, 2136, 1282, 2) and is consistent with io.fits
        self.shape = tuple([self.header['NAXIS{0}'.format(i + 1)]
                            for i in range(self.header['NAXIS'])][::-1])

        if len(self.shape) != 4:
            raise ValueError('Only 4D cube supported for now')

        self.strides = (
            self.shape[1] * self.shape[2] * self.shape[3] * self.itemsize,
            self.shape[2] * self.shape[3] * self.itemsize,
            self.shape[3] * self.itemsize,
            self.itemsize)

    def __getitem__(self, view):
        """``view`` is a tuple of slices."""

        # Not sure what @astrofrog means here:
        # TODO: normalize view to list of slices
        # TODO: take into account stride

        # TODO: what about BSCALE and BZERO?
        # TODO: hard-coded to 4D currently

        data_indiv = []
        start, end = view[3].start, view[3].stop
        i_stride_1 = start * self.strides[3]
        i_stride_2 = end * self.strides[3]

        for j in range(view[2].start, view[2].stop, view[2].step or 1):
            j_stride = j * self.strides[2]

            for k in range(view[1].start, view[1].stop, view[1].step or 1):
                k_stride = k * self.strides[1]

                for m in range(view[0].start, view[0].stop, view[0].step or 1):
                    m_stride = m * self.strides[0]

                    range_start = (self.data_start + i_stride_1 +
                                   j_stride + k_stride + m_stride)
                    range_end = (self.data_start + i_stride_2 +
                                 j_stride + k_stride + m_stride - 1)
                    content = self.get_data_single_range(
                        range_start, range_end)
                    data_indiv.append(content)

        data = b''.join(data_indiv)
        array = np.frombuffer(data, dtype=self.data_type)
        new_shape = (view[0].stop - view[0].start,
                     (view[1].stop - view[1].start) // (view[1].step or 1),
                     (view[2].stop - view[2].start) // (view[2].step or 1),
                     (view[3].stop - view[3].start) // (view[3].step or 1))

        return array.reshape(new_shape)

    def get_data_single_range(self, start, end):
        """Perform a GET request for single range.

        Parameters
        ----------
        start, end : int
            Starting and ending bytes to get.

        Returns
        -------
        content : str
            Bytes string from the request.

        """
        resp = self.s3_client.get_object(
            RequestPayer='requester', Bucket=self.s3_bucket, Key=self.s3_key,
            Range='bytes={}-{}'.format(start, end))
        return resp['Body'].read()


def subarray_by_bytes(fin, data_shape, data_type, view):
    """Offline version of the ``__get__`` method above for testing.

    .. note:: Make sure algorithm here is in-sync with ``__get__``.

    Parameters
    ----------
    fin : obj
        Binary file pointer to the FITS file, opened with
        ``open(filename, 'rb')``.

    data_shape : tuple
        Numpy array shape of the entire data array.

    data_type : obj
        Numpy data type for the array.

    view : tuple
        Tuple of slices.

    Returns
    -------
    subarray : ndarray
        Extracted subarray.

    """
    header_block_size = 2879
    data_start = 2 * header_block_size + 2
    itemsize = np.dtype(data_type).itemsize
    data_strides = (data_shape[1] * data_shape[2] * data_shape[3] * itemsize,
                    data_shape[2] * data_shape[3] * itemsize,
                    data_shape[3] * itemsize,
                    itemsize)
    print('strides:', data_strides)

    data_indiv = []
    start, end = view[3].start, view[3].stop
    i_stride_1 = start * data_strides[3]
    i_stride_2 = end * data_strides[3]

    for j in range(view[2].start, view[2].stop, view[2].step or 1):
        j_stride = j * data_strides[2]

        for k in range(view[1].start, view[1].stop, view[1].step or 1):
            k_stride = k * data_strides[1]

            for m in range(view[0].start, view[0].stop, view[0].step or 1):
                m_stride = m * data_strides[0]

                range_start = (data_start + i_stride_1 +
                               j_stride + k_stride + m_stride)
                range_end = (data_start + i_stride_2 +
                             j_stride + k_stride + m_stride - 1)

                print('Range={}-{}'.format(range_start, range_end))

                fin.seek(range_start)
                content = fin.read(range_end - range_start + 1)

                data_indiv.append(content)

    data = b''.join(data_indiv)
    array = np.frombuffer(data, dtype=data_type)
    print(array)

    new_shape = (view[0].stop - view[0].start,
                 (view[1].stop - view[1].start) // (view[1].step or 1),
                 (view[2].stop - view[2].start) // (view[2].step or 1),
                 (view[3].stop - view[3].start) // (view[3].step or 1))

    return array.reshape(new_shape)
