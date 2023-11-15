#!/usr/bin/env bash
#
# harvest_ce_wavss.sh
#
# Harvest the wavss data from all of the OOI Coastal Endurance moorings. Data
# sets include telemetered and recovered host. Data is downloaded from OOI Net
# and reworked to create a cleaner and more consistent set of files named and
# organized by the mooring, mooring sub-location, data delivery method and
# deployment.
#
# C. Wingard, 2019-12-16 -- Initial code

# set the base directory python command for all subsequent processing
. $(dirname $CONDA_EXE)/../etc/profile.d/conda.sh
conda activate ooi
PYTHON="python -m ooi_data_explorations.uncabled.process_wavss"

### CE02SHSM ###
BASE_FLAGS="-s CE02SHSM -n SBD12 -sn 05-WAVSSA000"
BASE_FILE="${HOME}/ooidata/m2m/ce02shsm/buoy/wavss/ce02shsm.buoy.wavss"
for i in $(seq -f "%02g" 1 17); do
  $PYTHON $BASE_FLAGS -mt telemetered -st wavss_a_dcl_statistics -dp $i -o "$BASE_FILE.deploy$i.telemetered.wavss_a_dcl_statistics.nc"
  $PYTHON $BASE_FLAGS -mt telemetered -st wavss_a_dcl_mean_directional -dp $i -o "$BASE_FILE.deploy$i.telemetered.wavss_a_dcl_mean_directional.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st wavss_a_dcl_statistics_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.wavss_a_dcl_statistics_recovered.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st wavss_a_dcl_mean_directional_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.wavss_a_dcl_mean_directional_recovered.nc"
done

### CE04OSSM ###
BASE_FLAGS="-s CE04OSSM -n SBD12 -sn 05-WAVSSA000"
BASE_FILE="${HOME}/ooidata/m2m/ce04ossm/buoy/wavss/ce04ossm.buoy.wavss"
for i in $(seq -f "%02g" 1 16); do
  $PYTHON $BASE_FLAGS -mt telemetered -st wavss_a_dcl_statistics -dp $i -o "$BASE_FILE.deploy$i.telemetered.wavss_a_dcl_statistics.nc"
  $PYTHON $BASE_FLAGS -mt telemetered -st wavss_a_dcl_mean_directional -dp $i -o "$BASE_FILE.deploy$i.telemetered.wavss_a_dcl_mean_directional.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st wavss_a_dcl_statistics_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.wavss_a_dcl_statistics_recovered.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st wavss_a_dcl_mean_directional_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.wavss_a_dcl_mean_directional_recovered.nc"
done

### CE07SHSM ###
BASE_FLAGS="-s CE07SHSM -n SBD12 -sn 05-WAVSSA000"
BASE_FILE="${HOME}/ooidata/m2m/ce07shsm/buoy/wavss/ce07shsm.buoy.wavss"
for i in $(seq -f "%02g" 1 17); do
  $PYTHON $BASE_FLAGS -mt telemetered -st wavss_a_dcl_statistics -dp $i -o "$BASE_FILE.deploy$i.telemetered.wavss_a_dcl_statistics.nc"
  $PYTHON $BASE_FLAGS -mt telemetered -st wavss_a_dcl_mean_directional -dp $i -o "$BASE_FILE.deploy$i.telemetered.wavss_a_dcl_mean_directional.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st wavss_a_dcl_statistics_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.wavss_a_dcl_statistics_recovered.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st wavss_a_dcl_mean_directional_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.wavss_a_dcl_mean_directional_recovered.nc"
done

### CE09OSSM ###
BASE_FLAGS="-s CE09OSSM -n SBD12 -sn 05-WAVSSA000"
BASE_FILE="${HOME}/ooidata/m2m/ce09ossm/buoy/wavss/ce09ossm.buoy.wavss"
for i in $(seq -f "%02g" 1 17); do
  $PYTHON $BASE_FLAGS -mt telemetered -st wavss_a_dcl_statistics -dp $i -o "$BASE_FILE.deploy$i.telemetered.wavss_a_dcl_statistics.nc"
  $PYTHON $BASE_FLAGS -mt telemetered -st wavss_a_dcl_mean_directional -dp $i -o "$BASE_FILE.deploy$i.telemetered.wavss_a_dcl_mean_directional.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st wavss_a_dcl_statistics_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.wavss_a_dcl_statistics_recovered.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st wavss_a_dcl_mean_directional_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.wavss_a_dcl_mean_directional_recovered.nc"
done
