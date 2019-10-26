import json

import boto3

import loadstraws as ls
import numpy as np


def lambda_handler(event, context):

    #Lambda Function to load a straw and take a nan-median of the image.
    
    #Inputs: TIC ID, Sector, camera, ccd, col, row
    bucket = "tess-straws"
    local_path = "/Users/smullally/TESS/tess-straws/"
    path = ""
    sector = 1
    camera = 1
    ccd = 1
    col = 221
    row = 250
    
    #cubeObj = ls.LoadTessCubeS3(bucket, path , "", "")
    cubeObj = ls.LoadTessCube(local_path)
    
    cube,cube_col,cube_row = cubeObj.get(camera, ccd, col, row, min_size_pix = 40)
    
    shape = np.shape(cube)
    
    print(shape)
    
    return {
        "statusCode": 200,
        "body": json.dumps({
                "shape": str(shape),
                "col": str(cube_col),
                "row": str(cube_row)
                })
        }

def test1():
    event = {}
    context = {}
    val = lambda_handler(event,context)
    print(val)
