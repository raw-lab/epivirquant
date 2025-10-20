#!/usr/bin/env bash

set -e

ENV_NAME=EpiVirQuant
echo "Creating conda environment: "$ENV_NAME

# initialize conda environment in bash script
eval "$(conda shell.bash hook)"

# create the metacerberus environment in conda
mamba create -y -n $ENV_NAME -c conda-forge joblib configargparse matplotlib numpy \
scipy scikit-image pypher stardist

conda activate $ENV_NAME

echo "Created conda environment: "$ENV_NAME
echo "Building packages"

pip install .

echo "Running EpiVirQuant with calibration image GSL-01"

epivirquant.py --dapi "data/GSL/tiff/dapi" --fitc "data/GSL/tiff/fitc" --calibration "data/GSL/tiff/dapi/GSL_+_blue_beads_1_(dapi).tiff"

echo "EpiVirQuant finished. To run EpiVirQuant on your own data change the --dapi, --fitc, and --calibration flag \
to where your data is. More information can be found by running epivirquant.py --help"