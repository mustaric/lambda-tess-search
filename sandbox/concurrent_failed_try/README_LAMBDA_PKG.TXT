Source: https://www.giladpeleg.com/blog/package-aws-lambda-python-code-with-docker/

Need to create these manually:

* Dockerfile -- Like the example, except delete "COPY lambdas ./"
* Pipfile -- Insert/edit Lambda dependencies here, like Numpy, Astropy, ...

Create ZIP file with these commands:

* sudo docker build . -t aws-lambda
* CONTAINER_ID=$(sudo docker create aws-lambda)
* sudo docker cp ${CONTAINER_ID}:/tmp/code.zip build.zip
* sudo docker rm ${CONTAINER_ID}

Clean the Docker when images and containers are not needed anymore:

* sudo docker images -a
* sudo docker rmi <ID> --force
* sudo docker container prune

Once you get "build.zip" in working directory:

* Copy it somewhere else. For example:
    * mkdir ~/foo
    * copy build.zip ~/foo
* Unzip it.
    * cd ~/foo
    * unzip build.zip
* Copy lambda module (e.g., "app_worker.py") to the top directory level.
    * cp /path/to/app_worker.py ~/foo

To re-zip it:

* cd ~/foo
* zip -r ../lambda_tess_worker.zip *

This "lambda_tess_worker.zip" file is what gets uploaded to AWS Lambda.
Similar workflow to deploy "app_main.py".
Some non-default AWS Lambda settings:

* Use Python 3.7.
* Give S3 bucket access. And also AWSLambdaRole for the Lambda calling Lambda.
* Increase timeout to 3s for "worker" and 15m for "main".
* Increase memory to 256 MB for "worker" and more (TBD) for "main".
* For Astropy caching to work, define a HOME env var and set it to /tmp .
