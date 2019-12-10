#!/usr/bin/env bash
#
# harvest_cp_phsen.sh
#
# Harvest the phsen data from all of the OOI Coastal Endurance moorings. Data
# sets include telemetered, recovered host and instrument data. Data is
# downloaded from OOI Net and reworked to create a cleaner and more consistent
# set of files named and organized by the mooring, mooring sub-location, data
# delivery method and deployment.
#
# C. Wingard, 2019-07-22 -- Initial code

# set the base directory python command for all subsequent processing
cd /home/ooiuser/code/ooi-data-explorations
PYTHON="python -m instruments.python.uncabled.phsen.request_phsen"

### CP01CNSM ###
BASE_FLAGS="-s CP01CNSM -n RID26 -sn 06-PHSEND000"
BASE_FILE="cp01cnsm/nsif/phsen/cp01cnsm.nsif.phsen"
for i in $(seq -f "%02g" 1 11); do
	$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.phsen_abcdef_dcl_instrument.nc"
	$PYTHON $BASE_FLAGS -mt recovered_host -st phsen_abcdef_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.phsen_abcdef_dcl_instrument_recovered.nc"
	$PYTHON $BASE_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_dcl_instrument -dp 12 -o "$BASE_FILE.deploy12.telemetered.phsen_abcdef_dcl_instrument.nc"

BASE_FLAGS="-s CP01CNSM -n MFD35 -sn 06-PHSEND000"
BASE_FILE="cp01cnsm/seafloor/phsen/cp01cnsm.seafloor.phsen"
for i in $(seq -f "%02g" 1 11); do
	$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.phsen_abcdef_dcl_instrument.nc"
	$PYTHON $BASE_FLAGS -mt recovered_host -st phsen_abcdef_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.phsen_abcdef_dcl_instrument_recovered.nc"
	$PYTHON $BASE_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_dcl_instrument -dp 12 -o "$BASE_FILE.deploy12.telemetered.phsen_abcdef_dcl_instrument.nc"

### CP03ISSM ###
BASE_FLAGS="-s CP03ISSM -n RID26 -sn 06-PHSEND000"
BASE_FILE="cp03issm/nsif/phsen/cp03issm.nsif.phsen"
for i in $(seq -f "%02g" 1 10); do
	$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.phsen_abcdef_dcl_instrument.nc"
	$PYTHON $BASE_FLAGS -mt recovered_host -st phsen_abcdef_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.phsen_abcdef_dcl_instrument_recovered.nc"
	$PYTHON $BASE_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_dcl_instrument -dp 11 -o "$BASE_FILE.deploy11.telemetered.phsen_abcdef_dcl_instrument.nc"

BASE_FLAGS="-s CP03ISSM -n MFD35 -sn 06-PHSEND000"
BASE_FILE="cp03issm/seafloor/phsen/cp03issm.seafloor.phsen"
for i in $(seq -f "%02g" 1 10); do
	$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.phsen_abcdef_dcl_instrument.nc"
	$PYTHON $BASE_FLAGS -mt recovered_host -st phsen_abcdef_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.phsen_abcdef_dcl_instrument_recovered.nc"
	$PYTHON $BASE_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_dcl_instrument -dp 11 -o "$BASE_FILE.deploy11.telemetered.phsen_abcdef_dcl_instrument.nc"

### CP04OSSM ###
BASE_FLAGS="-s CP04OSSM -n RID26 -sn 06-PHSEND000"
BASE_FILE="cp04ossm/nsif/phsen/cp04ossm.nsif.phsen"
for i in $(seq -f "%02g" 1 10); do
  $PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.phsen_abcdef_dcl_instrument.nc"
  $PYTHON $BASE_FLAGS -mt recovered_host -st phsen_abcdef_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.phsen_abcdef_dcl_instrument_recovered.nc"
  $PYTHON $BASE_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_dcl_instrument -dp 11 -o "$BASE_FILE.deploy11.telemetered.phsen_abcdef_dcl_instrument.nc"

BASE_FLAGS="-s CP04OSSM -n MFD35 -sn 06-PHSEND000"
BASE_FILE="cp04ossm/seafloor/phsen/cp04ossm.seafloor.phsen"
for i in $(seq -f "%02g" 1 10); do
    $PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_dcl_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.phsen_abcdef_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st phsen_abcdef_dcl_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.phsen_abcdef_dcl_instrument_recovered.nc"
    $PYTHON $BASE_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_dcl_instrument -dp 11 -o "$BASE_FILE.deploy11.telemetered.phsen_abcdef_dcl_instrument.nc"
