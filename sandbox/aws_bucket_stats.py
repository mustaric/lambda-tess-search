import boto3


def ffi_lc_cache_csv_stats(bucketname):
    """Find how many light curve data points generated by Lambda.
    This needs AWS credentials properly set up to work.

    Examples
    --------
    >>> d = ffi_lc_cache_csv_stats('mybucket')
    >>> d
    {'tic000025155310': {'s0001-4-1': 960,
      's0002-4-2': 979,
      's0003-4-2': 988,
      's0004-4-2': 870,
      's0005-4-3': 899,
      's0006-4-3': 814,
      's0007-4-3': 783,
      's0008-4-3': 700,
      's0009-4-4': 912,
      's0010-4-4': 915,
      's0011-4-4': 977,
      's0012-4-1': 960,
      's0013-4-1': 984},
     'tic000047525799': {'s0002-1-4': 960},
     'tic000070797900': {'s0002-1-2': 960},
     'tic000122220263': {'s0014-1-2': 884},
     'tic000144440290': {'s0002-2-4': 960},
     'tic000158276040': {'s0014-2-3': 913},
     'tic000201292545': {'s0001-2-2': 965},
     'tic000309792357': {'s0001-4-4': 952,
      's0002-4-1': 960,
      's0003-4-1': 988,
      's0005-4-2': 897,
      's0006-4-2': 814,
      's0007-4-2': 783,
      's0008-4-3': 700,
      's0009-4-3': 912,
      's0011-4-3': 977,
      's0012-4-4': 960,
      's0013-4-4': 984},
     'tic000388624270': {'s0008-3-4': 780, 's0009-3-3': 912},
     'tic000405454160': {'s0008-2-3': 700}}

    """
    s3 = boto3.resource('s3')
    bucket = s3.Bucket(bucketname)
    object_summary_iterator = bucket.objects.all()
    d = {}
    tot_uri = {}

    for obj in object_summary_iterator:
        key = obj.key

        # Skip items we do not care about.
        if not key.startswith('tic'):
            continue

        s = key.split('/')
        tic_id = s[0]
        sec_id = s[1]
        uri_list = f'tess-{sec_id}_s3_uris.txt'

        if sec_id not in tot_uri:
            s3_obj = s3.Object(bucketname, uri_list)
            content = s3_obj.get()['Body'].read()
            t = len(content.decode('utf8').split())
            tot_uri[sec_id] = t
        else:
            t = tot_uri[sec_id]

        if tic_id not in d:
            d[tic_id] = {}

        if sec_id in d[tic_id]:
            d[tic_id][sec_id][0] += 1
        else:
            d[tic_id][sec_id] = [1, t]

    return d
