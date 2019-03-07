#!/usr/bin/env bash

if [ -d /opt/conda/envs/rhaco ] && [ $BRANCH != 'master' ]; then
	echo "Using Cache";
	source activate rhaco;
else
	echo "Rebuilding Conda Env";
    rm -rf /opt/conda/envs/rhaco;
    conda create --name rhaco --clone base;
    source activate rhaco;
    pip install pytest-cov PyYAML codecov numba numpy-quaternion;
    pip install -e .;
fi
