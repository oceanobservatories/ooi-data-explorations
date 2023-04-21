#!/usr/bin/env bash
#
# harvest_ce_spkir.sh
#
# Harvest the SPKIR data from all of the OOI Coastal Endurance moorings. Data
# sets include telemetered and recovered host. Data is downloaded from OOI Net
# and reworked to create a cleaner and more consistent set of files named and
# organized by the mooring, mooring sub-location, data delivery method and
# deployment.
#
# C. Wingard, 2023-04-21 -- Initial code

# set the base directory python command for all subsequent processing
. $(dirname $CONDA_EXE)/../etc/profile.d/conda.sh
conda activate ooi
PYTHON="python -m ooi_data_explorations.uncabled.process_spkir"

### CE01ISSM ###
BASE_FLAGS="-s CE01ISSM -n RID16 -sn 08-SPKIRB000"
BASE_FILE="${HOME}/ooidata/m2m/ce01issm/nsif/spkir/ce01issm.nsif.spkir"
for i in $(seq -f "%02g" 1 18); do
    $PYTHON $BASE_FLAGS -mt telemetered -st spkir_abj_dcl_instrument -ba -dp $i -o "$BASE_FILE.deploy$i.telemetered.spkir_abj_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st spkir_abj_dcl_instrument_recovered -ba -dp $i -o "$BASE_FILE.deploy$i.recovered_host.spkir_abj_dcl_instrument_recovered.nc"
done

### CE02SHSM ###
BASE_FLAGS="-s CE02SHSM -n RID26 -sn 08-SPKIRB000"
BASE_FILE="${HOME}/ooidata/m2m/ce02shsm/nsif/spkir/ce02shsm.nsif.spkir"
for i in $(seq -f "%02g" 1 16); do
    $PYTHON $BASE_FLAGS -mt telemetered -st spkir_abj_dcl_instrument -ba -dp $i -o "$BASE_FILE.deploy$i.telemetered.spkir_abj_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st spkir_abj_dcl_instrument_recovered -ba -dp $i -o "$BASE_FILE.deploy$i.recovered_host.spkir_abj_dcl_instrument_recovered.nc"
done

### CE04OSSM ###
BASE_FLAGS="-s CE04OSSM -n RID26 -sn 08-SPKIRB000"
BASE_FILE="${HOME}/ooidata/m2m/ce04ossm/nsif/spkir/ce04ossm.nsif.spkir"
for i in $(seq -f "%02g" 1 15); do
    $PYTHON $BASE_FLAGS -mt telemetered -st spkir_abj_dcl_instrument -ba -dp $i -o "$BASE_FILE.deploy$i.telemetered.spkir_abj_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st spkir_abj_dcl_instrument_recovered -ba -dp $i -o "$BASE_FILE.deploy$i.recovered_host.spkir_abj_dcl_instrument_recovered.nc"
done

### CE06ISSM ###
BASE_FLAGS="-s CE06ISSM -n RID16 -sn 08-SPKIRB000"
BASE_FILE="${HOME}/ooidata/m2m/ce06issm/nsif/spkir/ce06issm.nsif.spkir"
for i in $(seq -f "%02g" 1 17); do
    $PYTHON $BASE_FLAGS -mt telemetered -st spkir_abj_dcl_instrument -ba -dp $i -o "$BASE_FILE.deploy$i.telemetered.spkir_abj_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st spkir_abj_dcl_instrument_recovered -ba -dp $i -o "$BASE_FILE.deploy$i.recovered_host.spkir_abj_dcl_instrument_recovered.nc"
done

### CE07SHSM ###
BASE_FLAGS="-s CE07SHSM -n RID26 -sn 08-SPKIRB000"
BASE_FILE="${HOME}/ooidata/m2m/ce07shsm/nsif/spkir/ce07shsm.nsif.spkir"
for i in $(seq -f "%02g" 1 16); do
    $PYTHON $BASE_FLAGS -mt telemetered -st spkir_abj_dcl_instrument -ba -dp $i -o "$BASE_FILE.deploy$i.telemetered.spkir_abj_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st spkir_abj_dcl_instrument_recovered -ba -dp $i -o "$BASE_FILE.deploy$i.recovered_host.spkir_abj_dcl_instrument_recovered.nc"
done

### CE09OSSM ###
BASE_FLAGS="-s CE09OSSM -n RID26 -sn 08-SPKIRB000"
BASE_FILE="${HOME}/ooidata/m2m/ce09ossm/nsif/spkir/ce09ossm.nsif.spkir"
for i in $(seq -f "%02g" 1 16); do
    $PYTHON $BASE_FLAGS -mt telemetered -st spkir_abj_dcl_instrument -ba -dp $i -o "$BASE_FILE.deploy$i.telemetered.spkir_abj_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st spkir_abj_dcl_instrument_recovered -ba -dp $i -o "$BASE_FILE.deploy$i.recovered_host.spkir_abj_dcl_instrument_recovered.nc"
done
