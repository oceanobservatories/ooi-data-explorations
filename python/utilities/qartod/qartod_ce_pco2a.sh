#!/usr/bin/env bash
#
# qartod_ce_pco2a.sh
#
# Collect the PCO2A data from the four OOI Coastal Endurance moorings to
# calculate QARTOD test ranges and generate the different lookup values and
# tables.
#
# C. Wingard, 2021-06-17 -- Initial code

# set the base directory python command for all subsequent processing
. $(dirname $CONDA_EXE)/../etc/profile.d/conda.sh
conda activate ooi
PYTHON="python -m ooi_data_explorations.qartod.endurance.qartod_ce_pco2a"

### CE02SHSM ###
$PYTHON -s CE02SHSM -n SBD12 -sn 04-PCO2AA000

### CE04OSSM ###
$PYTHON -s CE04OSSM -n SBD12 -sn 04-PCO2AA000

### CE07SHSM ###
$PYTHON -s CE07SHSM -n SBD12 -sn 04-PCO2AA000

### CE09OSSM ###
$PYTHON -s CE09OSSM -n SBD12 -sn 04-PCO2AA000
