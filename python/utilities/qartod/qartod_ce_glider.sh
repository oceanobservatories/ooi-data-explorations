#!/usr/bin/env bash
#
# qartod_ce_glider.sh
#
# Collect all CE glider data from the IOOS GliderDAC to calculate QARTOD
# test ranges and generate the gross range and climatology lookup tables
# for the CTD, DO, FLR and PAR sensors. For the ADCP, collect and process
# the data from the GC THREDDS catalog.
#
# C. Wingard, 2026-05-26 -- Initial code

. $(dirname $CONDA_EXE)/../etc/profile.d/conda.sh
conda activate ooi

# Use the GDAC for the CTD, DO, FLR and PAR sensors
python -m ooi_data_explorations.qartod.endurance.qartod_ce_glider_ctdgv -s CE05MOAS -n None -sn None
python -m ooi_data_explorations.qartod.endurance.qartod_ce_glider_dosta -s CE05MOAS -n None -sn None
python -m ooi_data_explorations.qartod.endurance.qartod_ce_glider_flort -s CE05MOAS -n None -sn None
python -m ooi_data_explorations.qartod.endurance.qartod_ce_glider_parad -s CE05MOAS -n None -sn None

## Use the GC THREDDS catalog for the ADCPA
#PYTHON="python -m ooi_data_explorations.qartod.endurance.qartod_ce_glider_adcpa"
#
#NODES=(G0871 G0917 G1012 G1134 G1135 G1136 G1137 G1153 \
#       GL247 GL311 GL312 GL319 GL320 GL326 GL327 \
#       GL381 GL382 GL383 GL384 GL386)
#
#mapfile -t SAMPLE < <(printf '%s\n' "${NODES[@]}" | shuf -n 5)
#
#for NODE in "${SAMPLE[@]}"; do
#    $PYTHON -s CE05MOAS -n $NODE -sn 03-ADCPAM000
#done
