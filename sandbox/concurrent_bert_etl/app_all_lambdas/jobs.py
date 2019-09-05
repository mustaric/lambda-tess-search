"""This module defines all the Lambda functions to be deployed to AWS.
The function name would be AWS Lambda function name as well. Everything to be
packaged to Lambda must be contained *inside* the function. Module level
imports are for local testing or packaging use only.

Use ``bert.constants.AWS_LAMBDA_FUNCTION`` to check if function is running
locally or in AWS.

"""
from bert import binding, constants, utils


@binding.follow('noop')
def bert_tess_fullframe_main_1():
    """Extract light curve data from TESS full frame images for a given
    TESS observation ID.

    This Lambda function loops through full frame images. For each full frame
    image, it invokes :func:`bert_tess_fullframe_worker`.

    """
    import os

    from astropy.coordinates import SkyCoord
    from astroquery.mast import Catalogs, Tesscut

    work_queue, done_queue, ologger = utils.comm_binders(
        bert_tess_fullframe_main_1)

    # TODO: Grab this from future higher level Lambda?
    #
    # NOTE: Hardcoding is necessary when DEBUG=false.
    # When DEBUG=true, can use manual test event but cannot use DynamoDB
    # streaming (i.e., worker will not be called automatically).
    if os.environ.get('DEBUG', 'true') == 'false':
        work_queue = [{
            "tic_id": "25155310",
            "radius": 2.5,
            "cutout_width": 30
        }]

    # https://exo.mast.stsci.edu/exomast_planet.html?planet=WASP126b
    # Example event:
    # {
    #   "tic_id": "25155310",
    #   "radius": 2.5,
    #   "cutout_width": 30
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
                'cutout_width': event['cutout_width']
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
    #   "cutout_width": 30
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
                'use_cache': 'false'})


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

    # Function to calculate circular mask.
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
        finally:
            # Clean up
            os.remove(outfilename)
