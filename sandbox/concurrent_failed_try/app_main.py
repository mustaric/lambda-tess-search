import json
# import os
from multiprocessing import Process, Pipe

import boto3
import numpy as np
from astropy import log
from astroquery.mast import Observations

# TODO: Remove this after https://github.com/astropy/astroquery/pull/1536
#       is merged.
log.setLevel('ERROR')

lambda_client = boto3.client('lambda')


def _pipe_worker(payload, conn):
    response = lambda_client.invoke(
        FunctionName='tess_fullframe_worker',
        InvocationType='Event',  # asynchronous
        Payload=payload)
    conn.send([response])
    conn.close()


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
    obs_id = event['id']  # TESS observation ID; Example: 'tess-s0001-1-1'

    # TODO: Calculate some of these from the 10th frame?
    # For now, also takes these and pass them onto worker:
    payload = {'xpos': event['xpos'],
               'ypos': event['ypos'],
               'radius': event['radius'],
               'bright_pixel_threshold': event['bright_pixel_threshold']}

    # Find full frame dataset for the observation ID.
    obs_table = Observations.query_criteria(obs_id=obs_id)
    products = Observations.get_product_list(obs_table)
    filtered = Observations.filter_products(
        products, productSubGroupDescription="FFIC", mrp_only=False)

    # Use AWS S3 bucket to pull data from.
    Observations.enable_cloud_dataset()  # TODO: verbose=False ?
    s3_urls = Observations.get_cloud_uris(filtered, include_bucket=False)

    # TODO: Timed out! Try https://docs.python.org/3/library/asyncio.html ?
    # TODO: Handle same Lambda call invoked multiple times by AWS?
    # Call tess_fullframe_worker AWS Lambda function in parallel
    # https://aws.amazon.com/blogs/compute/parallel-processing-in-python-with-aws-lambda/
    parent_connections = []
    processes = []
    data = []
    for url in s3_urls[:2]:  # TODO: Remove [:2] when done testing
        payload['key'] = url
        parent_conn, child_conn = Pipe()
        parent_connections.append(parent_conn)
        arg = json.dumps(payload)
        process = Process(target=_pipe_worker, args=(arg, child_conn))
        processes.append(process)

    for process in processes:
        process.start()

    for process in processes:
        process.join()

    for parent_connection in parent_connections:
        try:
            response = parent_connection.recv()[0]
        except EOFError:
            response = {}
        if 'body' not in response:  # Worker Lambda threw exception
            continue
        body = json.loads(response['body'])
        row = (body['midtime'], body['signal'], body['background'])
        if np.all(list(map(np.isfinite, row))):
            data.append(row)

    # TODO: Save data as table.
    # filename = f'/tmp/{obs_id}_lightcurve.csv'
    # with open(filename) as fout:
    #     for row in data:
    #         fout.write(f'{row[0]},{row[1]},{row[2]}{os.linesep}')

    # TODO: Upload table to S3 and then delete the table locally.
    # TODO: Return table S3 URL below.
    # TODO: Do we want to plot it and upload the plot too?
    #       If so, need to add matplotlib as dependency.

    return {"statusCode": 200,
            "body": json.dumps({'n_rows': len(data), 'data_url': 'TODO'})}
