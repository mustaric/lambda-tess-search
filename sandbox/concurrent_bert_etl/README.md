This directory contains code to be deployed to AWS with
[bert-etl](https://github.com/jbcurtin/bert-etl/).

AWS credentials must be set up in `~/.aws/credentials`.

Clone this repository and then `cd` into this directory.
Then, use `bert-etl` to deploy. Deploying this will use AWS Lambda and
DynamoDB, so charges might apply.

Example deployment command:

    bert-deploy.py -m app_all_lambdas.py
