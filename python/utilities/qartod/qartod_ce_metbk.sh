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
PYTHON="python -m ooi_data_explorations.qartod.endurance.qartod_ce_metbk"

### CE02SHSM ###
$PYTHON -s CE02SHSM -n SBD11 -sn 06-METBKA000 -co 2022-01-01T00:00:00

### CE04OSSM ###
$PYTHON -s CE04OSSM -n SBD11 -sn 06-METBKA000 -co 2022-01-01T00:00:00

### CE07SHSM ###
$PYTHON -s CE07SHSM -n SBD11 -sn 06-METBKA000 -co 2022-01-01T00:00:00

### CE09OSSM ###
$PYTHON -s CE09OSSM -n SBD11 -sn 06-METBKA000 -co 2022-01-01T00:00:00
