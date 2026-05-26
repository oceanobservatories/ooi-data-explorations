#!/usr/bin/env bash
#
# qartod_ce_fdchp.sh
#
# Collect the FDCHP data from the OOI Coastal Endurance Washington Shelf
# Surface Mooring to calculate QARTOD test ranges and generate the different
# lookup values and tables.
#
# C. Wingard, 2026-05-26 -- Initial code

# set the base directory python command for all subsequent processing
. $(dirname $CONDA_EXE)/../etc/profile.d/conda.sh
conda activate ooi
PYTHON="python -m ooi_data_explorations.qartod.endurance.qartod_ce_fdchp"

### CE02SHSM ###
$PYTHON -s CE02SHSM -n SBD12 -sn 08-FDCHPA000 -co 2026-01-01T00:00:00
