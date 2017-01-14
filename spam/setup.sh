#!/bin/sh

export SPAM_PATH=#Enter absolute path to spam folder
export SPAM_HOST=#Enter hostname of computer (command : hostname)
export PYTHON=#Enter pathname to Python installation (command : whereis python)
export PYTHONPATH=${SPAM_PATH}/python:${PYTHONPATH}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${SPAM_PATH}/lib
export PATH=${SPAM_PATH}/bin:${PATH}
