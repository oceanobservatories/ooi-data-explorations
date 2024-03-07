#!/usr/bin/env bash
#
# qartod_ce_metbk.sh
#
# Collect the METBK data from the four OOI Coastal Endurance moorings to
# calculate QARTOD test ranges and generate the different lookup values and
# tables.
#
# C. Wingard, 2022-09-08 -- Initial code

# set the base directory python command for all subsequent processing
. $(dirname $CONDA_EXE)/../etc/profile.d/conda.sh
conda activate ooi
PYTHON="python -m ooi_data_explorations.qartod.endurance.qartod_ce_vel3d"

### CE01ISSM ###
$PYTHON -s CE02SHSM -n MFD35 -sn 01-VEL3DD000 -co 2023-12-31T23:59:59.999

### CE06ISSM ###
$PYTHON -s CE04OSSM -n MFD35 -sn 01-VEL3DD000 -co 2023-12-31T23:59:59.999

### CE07SHSM ###
$PYTHON -s CE07SHSM -n MFD35 -sn 01-VEL3DD000 -co 2023-12-31T23:59:59.999

### CE09OSSM ###
$PYTHON -s CE09OSSM -n MFD35 -sn 01-VEL3DD000 -co 2023-12-31T23:59:59.999

### CE09OSPM ###
#$PYTHON -s CE09OSPM -n WFP01 -sn 01-VEL3DK000 -co 2023-12-31T23:59:59.999
