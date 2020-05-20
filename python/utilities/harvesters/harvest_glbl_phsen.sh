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
conda activate ooi
PYTHON="python -m ooi_data_explorations.uncabled.request_phsen"

### GA01SUMO ###
BASE_FLAGS="-s GA01SUMO -n RII11 -sn 02-PHSENE041"
BASE_FILE="${HOME}/ooidata/m2m/ga01sumo/imm/phsen/20m/ga01sumo.imm.phsen.20m"
for i in $(seq -f "%02g" 1 3); do
	$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_imodem_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.phsen_abcdef_imodem_instrument.nc"
	$PYTHON $BASE_FLAGS -mt recovered_host -st phsen_abcdef_imodem_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.phsen_abcdef_imodem_instrument_recovered.nc"
	$PYTHON $BASE_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done

BASE_FLAGS="-s GA01SUMO -n RII11 -sn 02-PHSENE042"
BASE_FILE="${HOME}/ooidata/m2m/ga01sumo/imm/phsen/100m/ga01sumo.imm.phsen.100m"
for i in $(seq -f "%02g" 1 3); do
	$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_imodem_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.phsen_abcdef_imodem_instrument.nc"
	$PYTHON $BASE_FLAGS -mt recovered_host -st phsen_abcdef_imodem_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.phsen_abcdef_imodem_instrument_recovered.nc"
	$PYTHON $BASE_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done

### GA03FLMA and GA03FLMB ###
FLMA_FLAGS="-s GA03FLMA -n RIS01 -sn 04-PHSENF000"
FLMA_FILE="ga03flma/nsif/phsen/ga03flma.nsif.phsen"
FLMB_FLAGS="-s GA03FLMB -n RIS01 -sn 04-PHSENF000"
FLMB_FILE="ga03flmb/nsif/phsen/ga03flmb.nsif.phsen"
for i in $(seq -f "%02g" 1 3); do
	$PYTHON $FLMA_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$FLMA_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
	$PYTHON $FLMB_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$FLMB_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done

### GI01SUMO ###
BASE_FLAGS="-s GI01SUMO -n RII11 -sn 02-PHSENE041"
BASE_FILE="${HOME}/ooidata/m2m/gi01sumo/imm/phsen/20m/gi01sumo.imm.phsen.20m"
for i in $(seq -f "%02g" 1 5); do
	$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_imodem_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.phsen_abcdef_imodem_instrument.nc"
	$PYTHON $BASE_FLAGS -mt recovered_host -st phsen_abcdef_imodem_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.phsen_abcdef_imodem_instrument_recovered.nc"
	$PYTHON $BASE_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_imodem_instrument -dp 6 -o "$BASE_FILE.deploy12.telemetered.phsen_abcdef_imodem_instrument.nc"

BASE_FLAGS="-s GI01SUMO -n RII11 -sn 02-PHSENE042"
BASE_FILE="${HOME}/ooidata/m2m/gi01sumo/imm/phsen/100m/gi01sumo.imm.phsen.100m"
for i in $(seq -f "%02g" 1 5); do
	$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_imodem_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.phsen_abcdef_imodem_instrument.nc"
	$PYTHON $BASE_FLAGS -mt recovered_host -st phsen_abcdef_imodem_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.phsen_abcdef_imodem_instrument_recovered.nc"
	$PYTHON $BASE_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_imodem_instrument -dp 6 -o "$BASE_FILE.deploy12.telemetered.phsen_abcdef_imodem_instrument.nc"

### GI03FLMA and GI03FLMB ###
FLMA_FLAGS="-s GI03FLMA -n RIS01 -sn 04-PHSENF000"
FLMA_FILE="gi03flma/nsif/phsen/gi03flma.nsif.phsen"
FLMB_FLAGS="-s GI03FLMB -n RIS01 -sn 04-PHSENF000"
FLMB_FILE="gi03flmb/nsif/phsen/gi03flmb.nsif.phsen"
for i in $(seq -f "%02g" 1 5); do
	$PYTHON $FLMA_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$FLMA_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
	$PYTHON $FLMB_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$FLMB_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done

### GP03FLMA and GP03FLMB ###
FLMA_FLAGS="-s GP03FLMA -n RIS01 -sn 04-PHSENF000"
FLMA_FILE="gp03flma/nsif/phsen/gp03flma.nsif.phsen"
FLMB_FLAGS="-s GP03FLMB -n RIS01 -sn 04-PHSENF000"
FLMB_FILE="gp03flmb/nsif/phsen/gp03flmb.nsif.phsen"
for i in $(seq -f "%02g" 1 6); do
	$PYTHON $FLMA_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$FLMA_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
	$PYTHON $FLMB_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$FLMB_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done

### GS01SUMO ###
BASE_FLAGS="-s GS01SUMO -n RII11 -sn 02-PHSENE041"
BASE_FILE="${HOME}/ooidata/m2m/gs01sumo/imm/phsen/20m/gs01sumo.imm.phsen.20m"
for i in $(seq -f "%02g" 1 3); do
	$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_imodem_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.phsen_abcdef_imodem_instrument.nc"
	$PYTHON $BASE_FLAGS -mt recovered_host -st phsen_abcdef_imodem_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.phsen_abcdef_imodem_instrument_recovered.nc"
	$PYTHON $BASE_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_imodem_instrument -dp 4 -o "$BASE_FILE.deploy12.telemetered.phsen_abcdef_imodem_instrument.nc"

BASE_FLAGS="-s GS01SUMO -n RII11 -sn 02-PHSENE042"
BASE_FILE="${HOME}/ooidata/m2m/gs01sumo/imm/phsen/100m/gs01sumo.imm.phsen.100m"
for i in $(seq -f "%02g" 1 3); do
	$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_imodem_instrument -dp $i -o "$BASE_FILE.deploy$i.telemetered.phsen_abcdef_imodem_instrument.nc"
	$PYTHON $BASE_FLAGS -mt recovered_host -st phsen_abcdef_imodem_instrument_recovered -dp $i -o "$BASE_FILE.deploy$i.recovered_host.phsen_abcdef_imodem_instrument_recovered.nc"
	$PYTHON $BASE_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st phsen_abcdef_imodem_instrument -dp 4 -o "$BASE_FILE.deploy12.telemetered.phsen_abcdef_imodem_instrument.nc"

### GS03FLMA and GS03FLMB ###
FLMA_FLAGS="-s GS03FLMA -n RIS01 -sn 04-PHSENF000"
FLMA_FILE="gs03flma/nsif/phsen/gs03flma.nsif.phsen"
FLMB_FLAGS="-s GS03FLMB -n RIS01 -sn 04-PHSENF000"
FLMB_FILE="gs03flmb/nsif/phsen/gs03flmb.nsif.phsen"
for i in $(seq -f "%02g" 1 3); do
	$PYTHON $FLMA_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$FLMA_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
	$PYTHON $FLMB_FLAGS -mt recovered_inst -st phsen_abcdef_instrument -dp $i -o "$FLMB_FILE.deploy$i.recovered_inst.phsen_abcdef_instrument.nc"
done
