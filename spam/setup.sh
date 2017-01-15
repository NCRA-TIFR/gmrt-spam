#!/bin/sh

export SPAM_PATH=/data2/shubhankar/gmrt-spam/spam
export SPAM_HOST=ISHWAR2
export PYTHON=/usr/bin/python2.7
export PYTHONPATH=${SPAM_PATH}/python:${PYTHONPATH}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${SPAM_PATH}/lib
export PATH=${SPAM_PATH}/bin:${PATH}
