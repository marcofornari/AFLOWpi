#!/bin/bash
set -e -u -x

function repair_wheel {
    wheel="$1"
    for whl in wheelhouse/*.whl; do
    auditwheel repair "$whl" -w /io/wheelhouse/
    done
    cp wheelhouse/*.whl /io/wheelhouse/
}

# Compile wheels
for PYBIN in /opt/python/*3*/bin; do
#    "${PYBIN}/pip" install -r /io/AFLOWpi/dev-requirements.txt
    "${PYBIN}/pip" install --upgrade pip
    "${PYBIN}/pip" wheel /io/ --no-deps -w wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    repair_wheel "$whl"
done

# Install packages and test
#for PYBIN in /opt/python/*/bin/; do

