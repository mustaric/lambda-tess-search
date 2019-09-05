This directory contains code to be deployed to AWS with
[bert-etl](https://github.com/jbcurtin/bert-etl/) (use dev version for now).

AWS credentials must be set up in `~/.aws/credentials`. Region
(N. Virginia, us-east-1) must be explicitly specified in the credentials file.

Need bert-etl.yaml file in the working directory where bert-deploy.py
is run.

Make sure to run this in a fresh Python 3.7 conda or virtual environment in
Linux because everything in site-packages will be packaged and uploaded to
Lambda. Need to install ipdb, boto3, requests, pyyaml, numpy, astropy, and
astroquery first.

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

For each Lambda function in jobs.py, all inputs come from `work_queue`, which
is populated by `done_queue` from the previous function. For the starter
function, its inputs come from `-i filename` in bert-deploy.py command.

To invoke the deployed Lambda without event arguments (the first one only
for now; you will not see logs on AWS site invoking it this way):

    bert-deploy.py -m app_all_lambdas -s aws-lambda -i

To un-deploy (delete deployed Lambda and DynamoDB on AWS; deleting DynamoDB can
take up to 10 minutes on AWS for large database):

    bert-deploy.py -m app_all_lambdas -s aws-lambda -u


To test locally
---------------

To run Lambda functions locally with Docker and redis without any AWS access
(except to download/upload data from S3), where `redis_docker`,
`docker-stop-all`, and `docker-remove-all` are custom bash commands to run
redis on Docker, to stop Docker processes, and to clean up Docker artifacts,
respectively:

    redis_docker

    bert-runner.py -m app_all_lambdas -f

    docker-stop-all

    docker-remove-all

    docker ps

To run Lambda functions locally but using AWS DynamoDB (charges may apply
for S3 downloads and DB access; do not use full loop!):

    BERT_QUEUE_TYPE=dynamodb bert-runner.py -m app_all_lambdas
