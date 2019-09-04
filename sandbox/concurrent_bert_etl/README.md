This directory contains code to be deployed to AWS with
[bert-etl](https://github.com/jbcurtin/bert-etl/).

AWS credentials must be set up in `~/.aws/credentials`. Region
(N. Virginia, us-east-1) must be explicitly specified in the credentials file.

Need bert-etl.yaml file in the working directory where bert-deploy.py
is run.

Make sure to run this in a fresh conda or virtual environment because
everything in site-packages will be packaged and uploaded to Lambda.
Need to install boto3, requests, pyyaml, numpy, astropy, astroquery first.

Clone this repository and then `cd` into this directory.
Then, use `bert-etl` to deploy. Deploying this will use AWS Lambda and
DynamoDB, so charges might apply.

Deployment command for local testing:

    bert-deploy.py -m app_all_lambdas -s aws-lambda -f

-f flushes the DynamoDB tables or the "redis" (local) tables. This option is
used when you only want to deploy DB and not Lambdas. This is a one-time
setup for local testing (see below).

Deployment command for production:

    bert-deploy.py -m app_all_lambdas -s aws-lambda


To test locally
---------------

To run Lambda functions locally but using AWS DynamoDB (charges may apply
for S3 downloads and DB access; do not use full loop!):

    BERT_QUEUE_TYPE=dynamodb bert-runner.py -m app_all_lambdas
