#!/usr/bin/env bash
#
# qartod_ce_parad.sh
#
# Collect the parad data from all of the OOI Coastal Endurance moorings and the
# CSPPs to calculate QARTOD test ranges and generate the different lookup
# values and tables.
#
# C. Wingard, 2022-06-17 -- Initial code

# set the base directory python command for all subsequent processing
. $(dirname $CONDA_EXE)/../etc/profile.d/conda.sh
conda activate ooi
PYTHON="python -m ooi_data_explorations.qartod.endurance.qartod_ce_parad"

### CE01ISSP ###
$PYTHON -s CE01ISSP -n SP001 -sn 10-PARADJ000 -co 2022-01-01T00:00:00

### CE02SHSP ###
$PYTHON -s CE02SHSP -n SP001 -sn 09-PARADJ000 -co 2022-01-01T00:00:00

### CE06ISSP ###
$PYTHON -s CE06ISSP -n SP001 -sn 10-PARADJ000 -co 2022-01-01T00:00:00

### CE07SHSP ###
$PYTHON -s CE07SHSP -n SP001 -sn 09-PARADJ000 -co 2022-01-01T00:00:00

### CE09OSPM ###
$PYTHON -s CE09OSPM -n WFP01 -sn 05-PARADK000 -co 2022-01-01T00:00:00
