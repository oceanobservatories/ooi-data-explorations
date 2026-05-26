#!/usr/bin/env bash
#
# qartod_ce_glider_adcpa.sh
#
# Collect the ADCPA data from the OOI Gold Copy THREDDS server for all CE
# gliders to calculate QARTOD test ranges and generate the gross range lookup
# table. Climatology tables are not generated for glider ADCPA data.
#
# A single representative reference designator is used; all CE glider ADCPA
# deployments are combined via the Gold Copy THREDDS catalog.
#
# C. Wingard, 2026-05-26 -- Initial code

. $(dirname $CONDA_EXE)/../etc/profile.d/conda.sh
conda activate ooi
PYTHON="python -m ooi_data_explorations.qartod.endurance.qartod_ce_glider_adcpa"

### CE05MOAS ###
$PYTHON -s CE05MOAS -n GL311 -sn 03-ADCPAM000
