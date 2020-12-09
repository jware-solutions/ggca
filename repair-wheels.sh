#!/bin/bash
export LD_LIBRARY_PATH=/usr/lib64/:$LD_LIBRARY_PATH

echo $LD_LIBRARY_PATH

for whl in target/wheels/*.whl; do
    auditwheel repair "$whl" -w ./target/wheels/wheelhouse
done