#!/usr/bin/env bash
#
# qartod_ce_presf.sh
#
# Collect the PRESF data from the four OOI Coastal Endurance moorings to
# calculate QARTOD test ranges and generate the different lookup values and
# tables.
#
# C. Wingard, 2021-06-17 -- Initial code

# set the base directory python command for all subsequent processing
. $(dirname $CONDA_EXE)/../etc/profile.d/conda.sh
conda activate ooi
PYTHON="python -m ooi_data_explorations.qartod.endurance.qartod_ce_presf"

### CE01ISSM ###
$PYTHON -s CE01ISSM -n MFD35 -sn 02-PRESFA000 -co 2021-01-01T00:00:00

### CE06ISSM ###
$PYTHON -s CE06ISSM -n MFD35 -sn 02-PRESFA000 -co 2021-01-01T00:00:00

### CE07SHSM ###
$PYTHON -s CE07SHSM -n MFD35 -sn 02-PRESFB000 -co 2021-01-01T00:00:00

### CE09OSSM ###
$PYTHON -s CE09OSSM -n MFD35 -sn 02-PRESFC000 -co 2021-01-01T00:00:00
