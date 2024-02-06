#!/bin/bash

source activate sister

pge_dir=$(cd "$(dirname "$0")" ; pwd -P)

echo "Creating runconfig"
python ${pge_dir}/generate_runconfig.py inputs.json

echo "Running L2a Spectral Resample PGE"
python ${pge_dir}/spectral_resample.py runconfig.json
