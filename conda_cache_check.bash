#!/usr/bin/env bash

if [ -d /opt/conda/envs/rhaco ] && [ $BRANCH != 'master' ]; then
    rm -rf /opt/conda/envs/rhaco
	conda env create -f environment.yml;
	source activate rhaco
else
	echo "Rebuilding Conda Env";
    rm -rf /opt/conda/envs/rhaco
	conda env create -f environment.yml;
        source activate rhaco
fi
