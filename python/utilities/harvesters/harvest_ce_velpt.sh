#!/usr/bin/env bash
#
# harvest_ce_velpt.sh
#
# Harvest the VELPT data from all of the OOI Coastal Endurance moorings. Data
# sets include telemetered and recovered host. Data is downloaded from OOI Net
# and reworked to create a cleaner and more consistent set of files named and
# organized by the mooring, mooring sub-location, data delivery method and
# deployment.
#
# C. Wingard, 2023-04-21 -- Initial code

# set the base directory python command for all subsequent processing
. $(dirname $CONDA_EXE)/../etc/profile.d/conda.sh
conda activate ooi
PYTHON="python -m ooi_data_explorations.uncabled.process_velpt"

### CE01ISSM ###
BASE_FLAGS="-s CE01ISSM -n SBD17 -sn 04-VELPTA000"
BASE_FILE="${HOME}/ooidata/m2m/ce01issm/buoy/velpt/ce01issm.buoy.velpt"
for i in $(seq -f "%02g" 1 18); do
    $PYTHON $BASE_FLAGS -mt telemetered -st velpt_ab_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.velpt_ab_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st velpt_ab_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_dcl_instrument_recovered.nc"
    $PYTHON $BASE_FLAGS -mt recovered_inst -st velpt_ab_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_instrument_recovered.nc"
done

BASE_FLAGS="-s CE01ISSM -n RID16 -sn 04-VELPTA000"
BASE_FILE="${HOME}/ooidata/m2m/ce01issm/nsif/velpt/ce01issm.nsif.velpt"
for i in $(seq -f "%02g" 1 18); do
    $PYTHON $BASE_FLAGS -mt telemetered -st velpt_ab_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.velpt_ab_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st velpt_ab_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_dcl_instrument_recovered.nc"
    $PYTHON $BASE_FLAGS -mt recovered_inst -st velpt_ab_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_instrument_recovered.nc"
done

### CE02SHSM ###
BASE_FLAGS="-s CE02SHSM -n SBD11 -sn 04-VELPTA000"
BASE_FILE="${HOME}/ooidata/m2m/ce02shsm/buoy/velpt/ce02shsm.buoy.velpt"
for i in $(seq -f "%02g" 1 16); do
    $PYTHON $BASE_FLAGS -mt telemetered -st velpt_ab_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.velpt_ab_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st velpt_ab_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_dcl_instrument_recovered.nc"
    $PYTHON $BASE_FLAGS -mt recovered_inst -st velpt_ab_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_instrument_recovered.nc"
done

BASE_FLAGS="-s CE02SHSM -n RID26 -sn 04-VELPTA000"
BASE_FILE="${HOME}/ooidata/m2m/ce02shsm/nsif/velpt/ce02shsm.nsif.velpt"
for i in $(seq -f "%02g" 1 16); do
    $PYTHON $BASE_FLAGS -mt telemetered -st velpt_ab_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.velpt_ab_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st velpt_ab_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_dcl_instrument_recovered.nc"
    $PYTHON $BASE_FLAGS -mt recovered_inst -st velpt_ab_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_instrument_recovered.nc"
done

### CE04OSSM ###
BASE_FLAGS="-s CE04OSSM -n SBD11 -sn 04-VELPTA000"
BASE_FILE="${HOME}/ooidata/m2m/ce04ossm/buoy/velpt/ce04ossm.buoy.velpt"
for i in $(seq -f "%02g" 1 15); do
    $PYTHON $BASE_FLAGS -mt telemetered -st velpt_ab_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.velpt_ab_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st velpt_ab_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_dcl_instrument_recovered.nc"
    $PYTHON $BASE_FLAGS -mt recovered_inst -st velpt_ab_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_instrument_recovered.nc"
done

BASE_FLAGS="-s CE04OSSM -n RID26 -sn 04-VELPTA000"
BASE_FILE="${HOME}/ooidata/m2m/ce04ossm/nsif/velpt/ce04ossm.nsif.velpt"
for i in $(seq -f "%02g" 1 15); do
    $PYTHON $BASE_FLAGS -mt telemetered -st velpt_ab_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.velpt_ab_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st velpt_ab_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_dcl_instrument_recovered.nc"
    $PYTHON $BASE_FLAGS -mt recovered_inst -st velpt_ab_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_instrument_recovered.nc"
done

### CE06ISSM ###
BASE_FLAGS="-s CE06ISSM -n SBD17 -sn 04-VELPTA000"
BASE_FILE="${HOME}/ooidata/m2m/ce06issm/buoy/velpt/ce06issm.buoy.velpt"
for i in $(seq -f "%02g" 1 17); do
    $PYTHON $BASE_FLAGS -mt telemetered -st velpt_ab_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.velpt_ab_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st velpt_ab_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_dcl_instrument_recovered.nc"
    $PYTHON $BASE_FLAGS -mt recovered_inst -st velpt_ab_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_instrument_recovered.nc"
done

BASE_FLAGS="-s CE06ISSM -n RID16 -sn 04-VELPTA000"
BASE_FILE="${HOME}/ooidata/m2m/ce06issm/nsif/velpt/ce06issm.nsif.velpt"
for i in $(seq -f "%02g" 1 17); do
    $PYTHON $BASE_FLAGS -mt telemetered -st velpt_ab_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.velpt_ab_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st velpt_ab_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_dcl_instrument_recovered.nc"
    $PYTHON $BASE_FLAGS -mt recovered_inst -st velpt_ab_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_instrument_recovered.nc"
done

### CE07SHSM ###
BASE_FLAGS="-s CE07SHSM -n SBD11 -sn 04-VELPTA000"
BASE_FILE="${HOME}/ooidata/m2m/ce07shsm/buoy/velpt/ce07shsm.buoy.velpt"
for i in $(seq -f "%02g" 1 16); do
    $PYTHON $BASE_FLAGS -mt telemetered -st velpt_ab_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.velpt_ab_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st velpt_ab_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_dcl_instrument_recovered.nc"
    $PYTHON $BASE_FLAGS -mt recovered_inst -st velpt_ab_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_instrument_recovered.nc"
done

BASE_FLAGS="-s CE07SHSM -n RID26 -sn 04-VELPTA000"
BASE_FILE="${HOME}/ooidata/m2m/ce07shsm/nsif/velpt/ce07shsm.nsif.velpt"
for i in $(seq -f "%02g" 1 16); do
    $PYTHON $BASE_FLAGS -mt telemetered -st velpt_ab_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.velpt_ab_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st velpt_ab_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_dcl_instrument_recovered.nc"
    $PYTHON $BASE_FLAGS -mt recovered_inst -st velpt_ab_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_instrument_recovered.nc"
done

### CE09OSSM ###
BASE_FLAGS="-s CE09OSSM -n SBD11 -sn 04-VELPTA000"
BASE_FILE="${HOME}/ooidata/m2m/ce09ossm/buoy/velpt/ce09ossm.buoy.velpt"
for i in $(seq -f "%02g" 1 16); do
    $PYTHON $BASE_FLAGS -mt telemetered -st velpt_ab_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.velpt_ab_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st velpt_ab_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_dcl_instrument_recovered.nc"
    $PYTHON $BASE_FLAGS -mt recovered_inst -st velpt_ab_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_instrument_recovered.nc"
done

BASE_FLAGS="-s CE09OSSM -n RID26 -sn 04-VELPTA000"
BASE_FILE="${HOME}/ooidata/m2m/ce09ossm/nsif/velpt/ce09ossm.nsif.velpt"
for i in $(seq -f "%02g" 1 16); do
    $PYTHON $BASE_FLAGS -mt telemetered -st velpt_ab_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.velpt_ab_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st velpt_ab_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_dcl_instrument_recovered.nc"
    $PYTHON $BASE_FLAGS -mt recovered_inst -st velpt_ab_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.velpt_ab_instrument_recovered.nc"
done
