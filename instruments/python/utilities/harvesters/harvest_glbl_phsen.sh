#!/usr/bin/env bash
#
# harvest_cp_phsen.sh
#
# Harvest the phsen data from all of the OOI Global moorings. Data
# sets include telemetered, recovered host and instrument data. Data is
# downloaded from OOI Net and reworked to create a cleaner and more consistent
# set of files named and organized by the mooring, mooring sub-location, data
# delivery method and deployment.
#
# C. Wingard, 2019-07-22 -- Initial code

# set the base directory python command for all subsequent processing
cd /home/ooiuser/code/ooi-data-explorations
PYTHON="python -m instruments.python.uncabled.request_phsen"

### GA01SUMO ###
BASE_FLAGS="-s GA01SUMO -n RII11 -sn 02-PHSENE041"
BASE_FILE="ga01sumo/imm/phsen/ga01sumo.imm.phsen.20m"
for i in $(seq -f "%02g" 1 3); do
	$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_imodem_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.phsen_abcdef_imodem_instrument.nc"
	$PYTHON $BASE_FLAGS -mt recovered_host -st phsen_abcdef_imodem_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.phsen_abcdef_imodem_instrument_recovered.nc"
	$PYTHON $BASE_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done

BASE_FLAGS="-s GA01SUMO -n RII11 -sn 02-PHSENE042"
BASE_FILE="ga01sumo/imm/phsen/ga01sumo.imm.phsen.100m"
for i in $(seq -f "%02g" 1 3); do
	$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_imodem_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.phsen_abcdef_imodem_instrument.nc"
	$PYTHON $BASE_FLAGS -mt recovered_host -st phsen_abcdef_imodem_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.phsen_abcdef_imodem_instrument_recovered.nc"
	$PYTHON $BASE_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done

### GI01SUMO ###
BASE_FLAGS="-s GI01SUMO -n RII11 -sn 02-PHSENE041"
BASE_FILE="gi01sumo/imm/phsen/gi01sumo.imm.phsen.20m"
for i in $(seq -f "%02g" 1 5); do
	$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_imodem_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.phsen_abcdef_imodem_instrument.nc"
	$PYTHON $BASE_FLAGS -mt recovered_host -st phsen_abcdef_imodem_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.phsen_abcdef_imodem_instrument_recovered.nc"
	$PYTHON $BASE_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_imodem_instrument -dp 12 -o "$BASE_FILE.deploy12.telemetered.phsen_abcdef_imodem_instrument.nc"

BASE_FLAGS="-s GI01SUMO -n RII11 -sn 02-PHSENE042"
BASE_FILE="gi01sumo/imm/phsen/gi01sumo.imm.phsen.100m"
for i in $(seq -f "%02g" 1 5); do
	$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_imodem_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.phsen_abcdef_imodem_instrument.nc"
	$PYTHON $BASE_FLAGS -mt recovered_host -st phsen_abcdef_imodem_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.phsen_abcdef_imodem_instrument_recovered.nc"
	$PYTHON $BASE_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_imodem_instrument -dp 12 -o "$BASE_FILE.deploy12.telemetered.phsen_abcdef_imodem_instrument.nc"

### GS01SUMO ###
BASE_FLAGS="-s GS01SUMO -n RII11 -sn 02-PHSENE041"
BASE_FILE="gs01sumo/imm/phsen/gs01sumo.imm.phsen.20m"
for i in $(seq -f "%02g" 1 3); do
	$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_imodem_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.phsen_abcdef_imodem_instrument.nc"
	$PYTHON $BASE_FLAGS -mt recovered_host -st phsen_abcdef_imodem_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.phsen_abcdef_imodem_instrument_recovered.nc"
	$PYTHON $BASE_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_imodem_instrument -dp 4 -o "$BASE_FILE.deploy12.telemetered.phsen_abcdef_imodem_instrument.nc"

BASE_FLAGS="-s GS01SUMO -n RII11 -sn 02-PHSENE042"
BASE_FILE="gs01sumo/imm/phsen/gs01sumo.imm.phsen.100m"
for i in $(seq -f "%02g" 1 5); do
	$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_imodem_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.phsen_abcdef_imodem_instrument.nc"
	$PYTHON $BASE_FLAGS -mt recovered_host -st phsen_abcdef_imodem_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.phsen_abcdef_imodem_instrument_recovered.nc"
	$PYTHON $BASE_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_imodem_instrument -dp 4 -o "$BASE_FILE.deploy12.telemetered.phsen_abcdef_imodem_instrument.nc"
