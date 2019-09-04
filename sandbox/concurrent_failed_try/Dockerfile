FROM amazonlinux:2

ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8
ENV PYTHONDONTWRITEBYTECODE 1

RUN yum -y update && yum -y install gcc python3-devel python3 zip git

RUN curl https://bootstrap.pypa.io/get-pip.py | python3 -
RUN pip install pipenv

WORKDIR /app
COPY Pipfile* ./

RUN pipenv install --deploy && cp -r $(pipenv --venv)/lib/python3.7/site-packages/. ./

RUN find . -type f -iname '*.pyc' -delete -o -type d -name __pycache__ -delete && rm Pipfile*

RUN zip -r /tmp/code.zip .
