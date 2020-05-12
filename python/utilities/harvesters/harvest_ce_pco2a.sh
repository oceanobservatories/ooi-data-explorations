#!/usr/bin/env bash
#
# harvest_cp_pco2a.sh
#
# Harvest the pco2a data from all of the OOI Coastal Endurance moorings. Data
# sets are limited to the telemetered. Data is downloaded from OOI Net and reworked
# to create a cleaner and more consistent set of files named and organized by the 
# mooring, mooring sub-location, data delivery method and deployment.
#
# C. Wingard, 2019-12-16 -- Initial code

# set the base directory python command for all subsequent processing
cd /home/ooiuser/code/ooi-data-explorations
OOI="/home/cwingard/anaconda3/envs/ooi/bin"
PYTHON="$OOI/python -m python.ooi_data_explorations.uncabled.request_pco2a"

### CE02SHSM ###
BASE_FLAGS="-s CE02SHSM -n SBD12 -sn 04-PCO2AA000"
BASE_FILE="ce02shsm/buoy/pco2a/ce02shsm.buoy.pco2a"
for i in $(seq -f "%02g" 1 10); do
  $PYTHON $BASE_FLAGS -mt telemetered -st pco2a_a_dcl_instrument_air -ba -dp $i -o "$BASE_FILE.deploy$i.telemetered.pco2a_a_dcl_instrument.nc"
  $PYTHON $BASE_FLAGS -mt telemetered -st pco2a_a_dcl_instrument_water -ba -dp $i -o "$BASE_FILE.deploy$i.telemetered.pco2a_a_dcl_instrument.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st pco2a_a_dcl_instrument_air_recovered -ba -dp $i -o "$BASE_FILE.deploy$i.recovered_host.pco2a_a_dcl_instrument_recovered.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st pco2a_a_dcl_instrument_water_recovered -ba -dp $i -o "$BASE_FILE.deploy$i.recovered_host.pco2a_a_dcl_instrument_recovered.nc"
done

### CE04OSSM ###
BASE_FLAGS="-s CE04OSSM -n SBD12 -sn 04-PCO2AA000"
BASE_FILE="ce04ossm/buoy/pco2a/ce04ossm.buoy.pco2a"
for i in $(seq -f "%02g" 1 9); do
  $PYTHON $BASE_FLAGS -mt telemetered -st pco2a_a_dcl_instrument_air -ba -dp $i -o "$BASE_FILE.deploy$i.telemetered.pco2a_a_dcl_instrument.nc"
  $PYTHON $BASE_FLAGS -mt telemetered -st pco2a_a_dcl_instrument_water -ba -dp $i -o "$BASE_FILE.deploy$i.telemetered.pco2a_a_dcl_instrument.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st pco2a_a_dcl_instrument_air_recovered -ba -dp $i -o "$BASE_FILE.deploy$i.recovered_host.pco2a_a_dcl_instrument_recovered.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st pco2a_a_dcl_instrument_water_recovered -ba -dp $i -o "$BASE_FILE.deploy$i.recovered_host.pco2a_a_dcl_instrument_recovered.nc"
done

### CE07SHSM ###
BASE_FLAGS="-s CE07SHSM -n SBD12 -sn 04-PCO2AA000"
BASE_FILE="ce07shsm/buoy/pco2a/ce07shsm.buoy.pco2a"
for i in $(seq -f "%02g" 1 10); do
  $PYTHON $BASE_FLAGS -mt telemetered -st pco2a_a_dcl_instrument_air -ba -dp $i -o "$BASE_FILE.deploy$i.telemetered.pco2a_a_dcl_instrument.nc"
  $PYTHON $BASE_FLAGS -mt telemetered -st pco2a_a_dcl_instrument_water -ba -dp $i -o "$BASE_FILE.deploy$i.telemetered.pco2a_a_dcl_instrument.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st pco2a_a_dcl_instrument_air_recovered -ba -dp $i -o "$BASE_FILE.deploy$i.recovered_host.pco2a_a_dcl_instrument_recovered.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st pco2a_a_dcl_instrument_water_recovered -ba -dp $i -o "$BASE_FILE.deploy$i.recovered_host.pco2a_a_dcl_instrument_recovered.nc"
done

### CE09OSSM ###
BASE_FLAGS="-s CE09OSSM -n SBD12 -sn 04-PCO2AA000"
BASE_FILE="ce09ossm/buoy/pco2a/ce09ossm.buoy.pco2a"
for i in $(seq -f "%02g" 1 10); do
  $PYTHON $BASE_FLAGS -mt telemetered -st pco2a_a_dcl_instrument_air -ba -dp $i -o "$BASE_FILE.deploy$i.telemetered.pco2a_a_dcl_instrument.nc"
  $PYTHON $BASE_FLAGS -mt telemetered -st pco2a_a_dcl_instrument_water -ba -dp $i -o "$BASE_FILE.deploy$i.telemetered.pco2a_a_dcl_instrument.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st pco2a_a_dcl_instrument_air_recovered -ba -dp $i -o "$BASE_FILE.deploy$i.recovered_host.pco2a_a_dcl_instrument_recovered.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st pco2a_a_dcl_instrument_water_recovered -ba -dp $i -o "$BASE_FILE.deploy$i.recovered_host.pco2a_a_dcl_instrument_recovered.nc"
done
