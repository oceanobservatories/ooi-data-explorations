#!/usr/bin/env bash
#
# harvest_cp_metbk.sh
#
# Harvest the metbk data from all of the OOI Coastal Endurance moorings. Data
# sets include telemetered and recovered host. Data is downloaded from OOI Net
# and reworked to create a cleaner and more consistent set of files named and
# organized by the mooring, mooring sub-location, data delivery method and
# deployment.
#
# C. Wingard, 2019-12-16 -- Initial code

# set the base directory python command for all subsequent processing
conda activate ooi
PYTHON="python -m ooi_data_explorations.uncabled.process_metbk"

### CE02SHSM ###
BASE_FLAGS="-s CE02SHSM -n SBD11 -sn 06-METBKA000"
BASE_FILE="${HOME}/ooidata/m2m/ce02shsm/buoy/metbk/ce02shsm.buoy.metbk"
for i in $(seq -f "%02g" 7 10); do
  $PYTHON $BASE_FLAGS -mt telemetered -st metbk_a_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.metbk_a_dcl_instrument.nc"
  $PYTHON $BASE_FLAGS -mt telemetered -st metbk_hourly -dp $i -o "$BASE_FILE.deploy$i.telemetered.hourly.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st metbk_a_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.metbk_a_dcl_instrument_recovered.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st metbk_hourly -dp $i -o "$BASE_FILE.deploy$i.recovered_host.hourly.nc"
done

### CE04OSSM ###
BASE_FLAGS="-s CE04OSSM -n SBD11 -sn 06-METBKA000"
BASE_FILE="${HOME}/ooidata/m2m/ce04ossm/buoy/metbk/ce04ossm.buoy.metbk"
for i in $(seq -f "%02g" 1 9); do
  $PYTHON $BASE_FLAGS -mt telemetered -st metbk_a_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.metbk_a_dcl_instrument.nc"
  $PYTHON $BASE_FLAGS -mt telemetered -st metbk_hourly -dp $i -o "$BASE_FILE.deploy$i.telemetered.hourly.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st metbk_a_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.metbk_a_dcl_instrument_recovered.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st metbk_hourly -dp $i -o "$BASE_FILE.deploy$i.recovered_host.hourly.nc"
done

### CE07SHSM ###
BASE_FLAGS="-s CE07SHSM -n SBD11 -sn 06-METBKA000"
BASE_FILE="${HOME}/ooidata/m2m/ce07shsm/buoy/metbk/ce07shsm.buoy.metbk"
for i in $(seq -f "%02g" 1 10); do
  $PYTHON $BASE_FLAGS -mt telemetered -st metbk_a_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.metbk_a_dcl_instrument.nc"
  $PYTHON $BASE_FLAGS -mt telemetered -st metbk_hourly -dp $i -o "$BASE_FILE.deploy$i.telemetered.hourly.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st metbk_a_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.metbk_a_dcl_instrument_recovered.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st metbk_hourly -dp $i -o "$BASE_FILE.deploy$i.recovered_host.hourly.nc"
done

### CE09OSSM ###
BASE_FLAGS="-s CE09OSSM -n SBD11 -sn 06-METBKA000"
BASE_FILE="${HOME}/ooidata/m2m/ce09ossm/buoy/metbk/ce09ossm.buoy.metbk"
for i in $(seq -f "%02g" 1 10); do
  $PYTHON $BASE_FLAGS -mt telemetered -st metbk_a_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.metbk_a_dcl_instrument.nc"
  $PYTHON $BASE_FLAGS -mt telemetered -st metbk_hourly -dp $i -o "$BASE_FILE.deploy$i.telemetered.hourly.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st metbk_a_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.metbk_a_dcl_instrument_recovered.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st metbk_hourly -dp $i -o "$BASE_FILE.deploy$i.recovered_host.hourly.nc"
done
