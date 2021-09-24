#!/bin/bash

#######################################################
## YOU NEED TO HAVE DOCKER INSTALLED AND RUN AS ROOT ##
#######################################################
rm -rf ./wheelhouse/*

export PLAT=manylinux2014_x86_64	
export DOCKER_IMAGE=quay.io/pypa/$PLAT

docker pull $DOCKER_IMAGE

docker run --rm -e PLAT=$PLAT -v `pwd`:/io $DOCKER_IMAGE $PRE_CMD /io/build-wheels.sh


export PLAT=manylinux_2_24_x86_64	
export DOCKER_IMAGE=quay.io/pypa/$PLAT

docker pull $DOCKER_IMAGE

docker run --rm -e PLAT=$PLAT -v `pwd`:/io $DOCKER_IMAGE $PRE_CMD /io/build-wheels.sh

#######################################################
## YOU NEED PYPI USERNAME AND PASSWORD FOR AFLOWPI!! ##
#######################################################

pip install twine

twine upload --repository-url https://test.pypi.org/legacy/ wheelhouse/*manylinux*.whl

###################################
## TO UPLOAD TO THE REAL PYPI    ##
## REMOVE THE "TEST." IN THE URL ##
###################################
