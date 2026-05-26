#!/usr/bin/env bash
#
# qartod_ce_flort.sh
#
# Collect the FLORT data from all of the OOI Coastal Endurance moorings and
# profilers to calculate QARTOD test ranges and generate the different lookup 
# values andtables.
#
# C. Wingard, 2022-02-16 -- Initial code

# set the base directory python command for all subsequent processing
. $(dirname $CONDA_EXE)/../etc/profile.d/conda.sh
conda activate ooi
PYTHON="python -m ooi_data_explorations.qartod.endurance.qartod_ce_flort"

### CE01ISSM ###
$PYTHON -s CE01ISSM -n SBD17 -sn 06-FLORTD000
$PYTHON -s CE01ISSM -n RID16 -sn 02-FLORTD000

### CE01ISSP ###
$PYTHON -s CE01ISSP -n SP001 -sn 07-FLORTJ000

### CE02SHSM ###
$PYTHON -s CE02SHSM -n RID27 -sn 02-FLORTD000

### CE02SHSP ###
$PYTHON -s CE02SHSP -n SP001 -sn 07-FLORTJ000

### CE04OSSM ###
$PYTHON -s CE04OSSM -n RID27 -sn 02-FLORTD000

### CE06ISSM ###
$PYTHON -s CE06ISSM -n SBD17 -sn 06-FLORTD000
$PYTHON -s CE06ISSM -n RID16 -sn 02-FLORTD000

### CE06ISSP ###
$PYTHON -s CE06ISSP -n SP001 -sn 07-FLORTJ000

### CE07SHSM ###
$PYTHON -s CE07SHSM -n RID27 -sn 02-FLORTD000

### CE07SHSP ###
$PYTHON -s CE07SHSP -n SP001 -sn 07-FLORTJ000

### CE09OSSM ###
$PYTHON -s CE09OSSM -n RID27 -sn 02-FLORTD000

### CE09OSPM ###
$PYTHON -s CE09OSPM -n WFP01 -sn 04-FLORTK000
