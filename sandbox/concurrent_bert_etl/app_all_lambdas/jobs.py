"""This module defines all the Lambda functions to be deployed to AWS.
The function name would be AWS Lambda function name as well. Everything to be
packaged to Lambda must be contained *inside* the function. Module level
imports are for local testing or packaging use only.

Use ``bert.constants.AWS_LAMBDA_FUNCTION`` to check if function is running
locally or in AWS.

When the ``DEBUG`` environment variable is set to ``true``, test event can
be run from Lambda dashboard, but it will not call the next Lambda because
it will not write to DynamoDB. When ``DEBUG`` is ``false``, the whole thing
will run but hardcoding the test case in the first function is necessary.

"""
from bert import binding, constants, utils


@binding.follow('noop')
def bert_tess_fullframe_main_0():
    """Farms out tasks for a list of TIC IDs."""

    import os

    import boto3

    work_queue, done_queue, ologger = utils.comm_binders(
        bert_tess_fullframe_main_0)

    s3 = boto3.resource('s3')
    homedir = os.environ.get('HOME')

    # TODO: How to run the full workflow without hardcoding the input here?
    # NOTE: See the note about DEBUG in module docstring.
    work_queue = [{"bucket": "ffi-lc-cache",
                   "key": "mullally_input_list_001.txt",
                   "use_cache": "true"}]

    # Example event:
    # {
    #   "bucket": "my_bucket",
    #   "key": "my_tic_id_list.csv",
    #   "use_cache": "true"
    # }
    for event in work_queue:
        bucketname = event['bucket']
        bucket = s3.Bucket(name=bucketname)
        key = event['key']
        basename = os.path.basename(key)
        filename = os.path.join(homedir, basename)

        try:
            ologger.info(f'Attempting to download {key} from {bucketname}')
            bucket.download_file(
                key, filename, ExtraArgs={"RequestPayer": "requester"})
        except Exception as exc:
            ologger.error(f'{basename}: {str(exc)}')
            continue
        else:
            ologger.info(f'Parsing {filename}')
            with open(filename) as fin:
                fin.readline()  # Skip header
                for line in fin:
                    row = line.split()[0].split(',')
                    tic_id = row[0]
                    ologger.info(f'Processing {tic_id}')
                    d = {'tic_id': tic_id,
                         'radius': float(row[1]),
                         'cutout_width': int(row[2]),
                         'use_cache': event['use_cache']}
                    done_queue.put(d)
        finally:
            # Clean up
            os.remove(filename)


@binding.follow(bert_tess_fullframe_main_0,
                pipeline_type=constants.PipelineType.CONCURRENT)
def bert_tess_fullframe_main_1():
    """Extract light curve data from TESS full frame images for a given
    TIC ID.

    This Lambda function loops through full frame images. For each full frame
    image, it invokes :func:`bert_tess_fullframe_worker`.

    """
    from astropy.coordinates import SkyCoord
    from astroquery.mast import Catalogs, Tesscut

    work_queue, done_queue, ologger = utils.comm_binders(
        bert_tess_fullframe_main_1)

    # https://exo.mast.stsci.edu/exomast_planet.html?planet=WASP126b
    # Example event:
    # {
    #   "tic_id": "25155310",
    #   "radius": 2.5,
    #   "cutout_width": 30,
    #   "use_cache": "true"
    # }
    for event in work_queue:
        tic_id = event['tic_id']
        info = Catalogs.query_criteria(catalog='Tic', ID=tic_id)
        ra = info['ra'][0]  # deg
        dec = info['dec'][0]  # deg
        coo = SkyCoord(ra, dec, unit='deg')

        # TODO: Allow auto determination of radius from mag if not given.
        # mag = info['Tmag'][0]

        sector_resp = Tesscut.get_sectors(coordinates=coo)
        ologger.info(f'Found {len(sector_resp)} sectors for {tic_id}')

        # Process all the matched sectors.
        for sec_id in sector_resp['sectorName']:
            done_queue.put({
                'tic_id': tic_id,
                'sec_id': sec_id,
                'ra': ra,
                'dec': dec,
                'radius': event['radius'],
                'cutout_width': event['cutout_width'],
                'use_cache': event['use_cache']
            })


@binding.follow(bert_tess_fullframe_main_1,
                pipeline_type=constants.PipelineType.CONCURRENT)
def bert_tess_fullframe_main_2():
    """Continuation of main function to run it across different sectors."""

    import os
    import time

    import boto3
    from astropy.io import fits
    from astropy.wcs import WCS
    from astroquery.mast import Observations

    s3 = boto3.resource('s3')
    bucket = s3.Bucket(name=os.environ.get('AWSBUCKETNAME'))
    outbucket = s3.Bucket(name=os.environ.get('CACHEBUCKETNAME'))
    homedir = os.environ.get('HOME')

    work_queue, done_queue, ologger = utils.comm_binders(
        bert_tess_fullframe_main_2)

    # Example event:
    # {
    #   "tic_id": "25155310",
    #   "sec_id": "tess-s0001-4-1",
    #   "ra": 63.3739396231274,
    #   "dec": -69.226822697583,
    #   "radius": 2.5,
    #   "cutout_width": 30,
    #   "use_cache": "true"
    # }
    #
    # work_queue populated by calling Lambda
    for event in work_queue:
        tic_id = event['tic_id']
        sec_id = event['sec_id']

        basename = f'{sec_id}_s3_uris.txt'  # noqa
        filename = os.path.join(homedir, basename)

        try:
            # Check if URI list already cached.
            # According to MAST, there is no need to invalidate cache here.
            ologger.info(f'Attempting to download {basename} from S3')
            outbucket.download_file(
                basename, filename,
                ExtraArgs={"RequestPayer": "requester"})
        except Exception:
            # Find full frame dataset for the observation ID.
            ologger.info('Started quering Observations...')
            obs_table = Observations.query_criteria(obs_id=sec_id)
            products = Observations.get_product_list(obs_table)
            filtered = Observations.filter_products(
                products, productSubGroupDescription="FFIC",
                mrp_only=False)

            # Use AWS S3 bucket to pull data from.
            Observations.enable_cloud_dataset(verbose=False)
            ologger.info('Started obtaining cloud URIs...')
            t_start = time.time()
            s3_urls = Observations.get_cloud_uris(
                filtered, include_bucket=False)
            t_end = time.time()
            ologger.info(f'Got {len(s3_urls)} URIs in {t_end - t_start} s')

            # Upload URI list to cache.
            with open(filename, 'w') as fout:
                for url in s3_urls:
                    fout.write(url + os.linesep)
            try:
                outbucket.upload_file(
                    filename, basename,
                    ExtraArgs={"RequestPayer": "requester"})
            except Exception as exc:
                ologger.error(str(exc))
            else:
                ologger.info(f'Uploaded {basename} to S3')
        else:
            # Use cache if it exists.
            with open(filename, 'r') as fin:
                s3_urls = [url.strip() for url in fin.readlines()]
            ologger.info(f'Read {len(s3_urls)} URIs from {basename}')
        finally:
            # Clean up
            if os.path.exists(filename):
                os.remove(filename)

        ra = float(event['ra'])
        dec = float(event['dec'])

        # Find pixel coordinates from sky from first frame header.
        key = s3_urls[0]
        basename = key.split('/')[-1]
        filename = os.path.join(homedir, basename)
        ologger.info(f'Resolving WCS from {key}')
        bucket.download_file(
            key, filename, ExtraArgs={"RequestPayer": "requester"})
        hdr = fits.getheader(filename, ext=1)
        w = WCS(hdr)
        pix = w.all_world2pix(ra, dec, 0)
        xpos = round(float(pix[0]))  # float needed to get rid of 0-D array
        ypos = round(float(pix[1]))

        # Clean up
        os.remove(filename)

        # The star needs to be at least 2*radii pixels away in both X and Y.
        radius = float(event['radius'])
        edge_r = 2 * radius
        naxis1, naxis2 = w.pixel_shape  # X Y
        if (xpos < edge_r or xpos >= (naxis1 - edge_r) or
                ypos < edge_r or ypos >= (naxis2 - edge_r)):
            ologger.error(
                f'TIC f{tic_id} in {sec_id}: X={xpos},Y={ypos} not at least '
                f'{edge_r} pixels away from the edge, skipping...')
            continue

        # Pass data into the next AWS Lambda function.
        ologger.info(f'TIC f{tic_id} in {sec_id}: Started processing '
                     'full frame URIs...')
        for url in s3_urls:
            done_queue.put({
                'key': url,
                'tic_id': tic_id,
                'ra': ra,
                'dec': dec,
                'xpos': xpos,
                'ypos': ypos,
                'radius': radius,
                'cutout_width': event['cutout_width'],
                'use_cache': event['use_cache']})


@binding.follow(bert_tess_fullframe_main_2,
                pipeline_type=constants.PipelineType.CONCURRENT)
def bert_tess_fullframe_worker_1():
    """Extract light curve data from one TESS full frame image.

    This function is meant to be called from
    :func:`bert_tess_fullframe_main_2` concurrently.

    """
    import os

    import boto3
    import numpy as np
    from astropy.io import fits

    np.seterr(all='ignore')

    # Function to calculate circular and square masks.
    # https://stackoverflow.com/questions/49330080/numpy-2d-array-selecting-indices-in-a-circle
    def get_masks(arr_shape, cx, cy, r, w):
        y = np.arange(arr_shape[0])
        x = np.arange(arr_shape[1])
        circ = ((x[np.newaxis, :] - cx) ** 2 + (y[:, np.newaxis] - cy) ** 2) < (r ** 2)  # noqa
        hw = w // 2
        sqr = (np.abs(x[np.newaxis, :] - cx) <= hw) & (np.abs(y[:, np.newaxis] - cy) <= hw)  # noqa
        return circ, sqr

    s3 = boto3.resource('s3')
    bucket = s3.Bucket(name=os.environ.get('AWSBUCKETNAME'))
    outbucket_name = os.environ.get('CACHEBUCKETNAME')
    outbucket = s3.Bucket(name=outbucket_name)
    homedir = os.environ.get('HOME')

    work_queue, done_queue, ologger = utils.comm_binders(
        bert_tess_fullframe_worker_1)

    # Example event:
    # {
    #   "key": "tess/public/ffi/s0001/2018/206/4-1/tess2018206192942-s0001-4-1-0120-s_ffic.fits",  # noqa
    #   "tic_id": "25155310",
    #   "ra": 63.3739396231274,
    #   "dec": -69.226822697583,
    #   "xpos": 1623,
    #   "ypos": 392,
    #   "radius": 2.5,
    #   "cutout_width": 30,
    #   "use_cache": "true"
    # }
    #
    # work_queue populated by calling Lambda
    for event in work_queue:
        key = event['key']
        basename = key.split('/')[-1]
        filename = os.path.join(homedir, basename)
        sec_id = '-'.join(basename.split('-')[1:4])
        tic_id = event['tic_id']
        radius = float(event['radius'])
        cutout_width = int(event['cutout_width'])

        # Use cached CSV generated by previous run and skip recalculations.
        use_cache = event['use_cache'] == 'true'

        # Output CSV containing one data point of the light curve.
        # Example s3key: ticXXXXX/s0001-4-1/r2.5/w30/blahblah.csv
        outbasename = basename.replace('.fits', '.csv')
        outfilename = os.path.join(homedir, outbasename)
        s3key = f'tic{tic_id:0>12}/{sec_id}/r{radius}/w{cutout_width}/{outbasename}'  # noqa

        # If this output exists and user wants to use the cache, there is
        # nothing to do.
        if use_cache:
            try:
                s3.Object(outbucket_name, s3key).load()
            except Exception:  # Does not exist
                pass
            else:  # It exists; nothing to do
                ologger.info(f'{s3key} exists, skipping...')
                continue

        # TODO: Find a way to not do this when testing locally to save money?
        ologger.info(f'Downloading {key}')
        try:
            bucket.download_file(
                key, filename, ExtraArgs={"RequestPayer": "requester"})
        except Exception as exc:  # Download failed
            ologger.error(f'{key}: {str(exc)}')
            if os.path.exists(filename):  # Clean up
                os.remove(filename)
            continue

        ologger.info(f'Calculating data point for TIC {tic_id} in {basename}')
        xpos = int(event['xpos'])
        ypos = int(event['ypos'])
        ra = float(event['ra'])
        dec = float(event['dec'])

        with fits.open(filename) as pf:
            # EXT 0: Primary header
            # Observation time in BJD.
            obstime = 0.5 * (pf[0].header['TSTOP'] + pf[0].header['TSTART'])

            # EXT 1: cal
            dqflag = pf[1].header['DQUALITY']

            circ_mask, sqr_mask = get_masks(
                pf[1].data.shape, xpos, ypos, radius, cutout_width)

            # Signal is inside a circle.
            sig_arr = pf[1].data[circ_mask]

            # Sky background is outside circle but inside box.
            bg_arr = pf[1].data[~circ_mask & sqr_mask]
            sky = np.nanmedian(bg_arr)

            # Simple aperture photometry.
            signal = np.nansum(sig_arr - sky)

        # Clean up
        os.remove(filename)

        # Write data point to CSV
        ologger.info(f'TIC {tic_id}: Writing {outfilename}')
        with open(outfilename, 'w') as fout:
            fout.write(f'{obstime},{signal},{sky},{dqflag},'
                       f'{xpos},{ypos},{ra},{dec}')

        # Upload CSV to S3 bucket
        try:
            outbucket.upload_file(outfilename, s3key,
                                  ExtraArgs={"RequestPayer": "requester"})
        except Exception as exc:
            ologger.error(str(exc))
        else:
            ologger.info(f'Uploaded {s3key} to S3')
            done_queue.put({
                'tic_id': tic_id,
                'sec_id': sec_id,
                'radius': radius,
                'cutout_width': cutout_width,
                'use_cache': event['use_cache']})
        finally:
            # Clean up
            os.remove(outfilename)


@binding.follow(bert_tess_fullframe_worker_1,
                pipeline_type=constants.PipelineType.BOTTLE)
def bert_tess_fullframe_worker_2():
    """Collect light curve data from individual full frame images.
    Use the data to build a light curve file.
    Then, upload the file to S3 bucket.

    .. note:: Do not deploy from ``concurrent_bert_etl_worker2``
              folder anymore.

    TESS FFI Light Curve Format documented at
    https://archive.stsci.edu/missions/tess/doc/EXP-TESS-ARC-ICD-TM-0014.pdf#page=32

    """
    import os
    from datetime import datetime

    import boto3
    from astropy.table import Table

    work_queue, done_queue, ologger = utils.comm_binders(
        bert_tess_fullframe_worker_2)

    s3 = boto3.resource('s3')
    inbucket = s3.Bucket(name=os.environ.get('CACHEBUCKETNAME'))
    bucket_name = os.environ.get('AWSBUCKETNAME')
    bucket = s3.Bucket(name=bucket_name)
    homedir = os.environ.get('HOME')

    # Example event:
    # {
    #   "tic_id": "25155310",
    #   "sec_id": "s0001-4-1"
    #   "radius": 2.5,
    #   "cutout_width": 30,
    #   "use_cache": "true"
    # }
    #
    # work_queue populated by calling Lambda
    for event in work_queue:
        tic_id = event['tic_id']
        sec_id = event['sec_id']
        radius = float(event['radius'])
        cutout_width = int(event['cutout_width'])

        in_pfx = f'tic{tic_id:0>12}/{sec_id}/r{radius}/w{cutout_width}'
        basename = f'tic{tic_id:0>12}_{sec_id}_lcc.fits'
        s3key = f'tic{tic_id:0>12}/{basename}'
        outfilename = os.path.join(homedir, basename)

        # Use cached LC generated by previous run and skip recalculations.
        # Skipping also means BLS Lambda listening for S3 upload will not run.
        use_cache = event['use_cache'] == 'true'

        # If this output exists and user wants to use the cache, there is
        # nothing to do.
        if use_cache:
            try:
                s3.Object(bucket_name, s3key).load()
            except Exception:  # Does not exist
                pass
            else:  # It exists; nothing to do
                ologger.info(f'{s3key} exists, skipping...')
                continue

        sec_id_split = sec_id.split('-')
        sector = int(sec_id_split[0][1:])
        camera = int(sec_id_split[1])
        ccd = int(sec_id_split[2])

        # Table header
        lc_meta = {
            'TELESCOP': 'TESS',
            'CAMERA': camera,
            'SECTOR': sector,
            'CCD': ccd,
            'OBJECT': f'TIC {tic_id}',
            'RADESYS': 'ICRS',
            'AP_RAD': radius,
            'SKYWIDTH': cutout_width,
            'DATE': datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')}

        # f4 = np.float32, f8 = np.float64, i4 = np.int32
        lc_tab = Table(names=('TIME', 'SAP_FLUX', 'SAP_BKG', 'QUALITY'),
                       dtype=('f8', 'f4', 'f4', 'i4'),
                       meta=lc_meta)

        # Grab all the light curve data points and piece them together.
        for obj in inbucket.objects.filter(
                Prefix=in_pfx, RequestPayer='requester'):
            filename = os.path.join(homedir, os.path.basename(obj.key))
            inbucket.download_file(
                obj.key, filename, ExtraArgs={"RequestPayer": "requester"})

            with open(filename, 'r') as fin:
                row = fin.read().split(',')

            # Clean up
            os.remove(filename)

            midtime = float(row[0])
            signal = float(row[1])
            background = float(row[2])
            dqflag = int(row[3])
            xpos = int(row[4])
            ypos = int(row[5])
            ra = float(row[6])
            dec = float(row[7])

            lc_tab.add_row((midtime, signal, background, dqflag))

        # Sort table by observation time.
        lc_tab.sort('TIME')

        # More metadata
        lc_tab.meta.update({
            'RA_OBJ': ra,
            'DEC_OBJ': dec,
            'APCEN_X': xpos,
            'APCEN_Y': ypos})

        # Write locally to FITS table.
        # Table data and metadata will go to EXT 1.
        lc_tab.write(outfilename, format='fits')
        ologger.info(f'Light curve [{outfilename}] with {len(lc_tab)} pts')

        # Upload to S3 bucket.
        try:
            bucket.upload_file(
                outfilename, s3key, ExtraArgs={"RequestPayer": "requester"})
        except Exception as exc:
            ologger.error(str(exc))
        else:
            ologger.info(f'Uploaded {s3key} to S3')
        finally:
            # Clean up
            os.remove(outfilename)
