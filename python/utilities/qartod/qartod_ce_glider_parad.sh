#!/usr/bin/env bash
#
# qartod_ce_glider_parad.sh
#
# Collect all CE glider PAR data from the IOOS GliderDAC to calculate QARTOD
# test ranges and generate the gross range and climatology lookup tables.
#
# A single representative reference designator is used; the GliderDAC search
# covers the full CE array (extent=100 nm), so all glider deployments are
# included regardless of which GL node is specified.
#
# C. Wingard, 2026-05-26 -- Initial code

. $(dirname $CONDA_EXE)/../etc/profile.d/conda.sh
conda activate ooi
PYTHON="python -m ooi_data_explorations.qartod.endurance.qartod_ce_glider_parad"

### CE05MOAS ###
$PYTHON -s CE05MOAS -n GL311 -sn 01-PARADM000
