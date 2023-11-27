#!/usr/bin/env bash
#
# generate_fcoeff.sh
#
# Use the Log 9 formatted bulk wave data to generate the frequency
# coefficients for a mooring and deployment of interest.
#
# C. Wingard, 2023-11-27 -- Initial code

# make sure the user has provided a Log 9 file to process
if [ $# -eq 0 ]; then
    echo "No Log 9 file provided"
    exit 1
fi
LOG9_FILE="$1"

# set the base directory python command for all subsequent processing
. $(dirname "$CONDA_EXE")/../etc/profile.d/conda.sh
conda activate ooi

# process the Log 9 data, creating the frequency coefficients
python -m ooi_data_explorations.uncabled.utilities.generate_fcoeff -l "$LOG9_FILE"