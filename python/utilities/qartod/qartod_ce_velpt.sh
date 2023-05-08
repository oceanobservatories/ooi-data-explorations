#!/usr/bin/env bash
#
# qartod_ce_flort.sh
#
# Collect the VELPT data from all of the OOI Coastal Endurance moorings and
# profilers to calculate QARTOD test ranges and generate the different lookup 
# values andtables.
#
# C. Wingard, 2023-02-16 -- Initial code

# set the base directory python command for all subsequent processing
. $(dirname $CONDA_EXE)/../etc/profile.d/conda.sh
conda activate ooi
PYTHON="python -m ooi_data_explorations.qartod.endurance.qartod_ce_velpt"

### CE01ISSM ###
$PYTHON -s CE01ISSM -n SBD17 -sn 04-VELPTA000 -co 2023-01-01T00:00:00
$PYTHON -s CE01ISSM -n RID16 -sn 04-VELPTA000 -co 2023-01-01T00:00:00

### CE01ISSP ###
$PYTHON -s CE01ISSP -n SP001 -sn 05-VELPTJ000 -co 2023-01-01T00:00:00

### CE02SHSM ###
$PYTHON -s CE02SHSM -n SBD11 -sn 04-VELPTA000 -co 2023-01-01T00:00:00
$PYTHON -s CE02SHSM -n RID26 -sn 04-VELPTA000 -co 2023-01-01T00:00:00

### CE02SHSP ###
$PYTHON -s CE02SHSP -n SP001 -sn 02-VELPTJ000 -co 2023-01-01T00:00:00

### CE04OSSM ###
$PYTHON -s CE04OSSM -n SBD11 -sn 04-VELPTA000 -co 2023-01-01T00:00:00
$PYTHON -s CE04OSSM -n RID26 -sn 04-VELPTA000 -co 2023-01-01T00:00:00

### CE06ISSM ###
$PYTHON -s CE06ISSM -n SBD17 -sn 04-VELPTA000 -co 2023-01-01T00:00:00
$PYTHON -s CE06ISSM -n RID16 -sn 04-VELPTA000 -co 2023-01-01T00:00:00

### CE06ISSP ###
$PYTHON -s CE06ISSP -n SP001 -sn 05-VELPTJ000 -co 2023-01-01T00:00:00

### CE07SHSM ###
$PYTHON -s CE07SHSM -n SBD11 -sn 04-VELPTA000 -co 2023-01-01T00:00:00
$PYTHON -s CE07SHSM -n RID26 -sn 04-VELPTA000 -co 2023-01-01T00:00:00

### CE07SHSP ###
$PYTHON -s CE07SHSP -n SP001 -sn 02-VELPTJ000 -co 2023-01-01T00:00:00

### CE09OSSM ###
$PYTHON -s CE09OSSM -n SBD11 -sn 04-VELPTA000 -co 2023-01-01T00:00:00
$PYTHON -s CE09OSSM -n RID26 -sn 04-VELPTA000 -co 2023-01-01T00:00:00
