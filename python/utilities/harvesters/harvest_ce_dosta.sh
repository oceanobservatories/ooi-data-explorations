#!/usr/bin/env bash
#
# harvest_ce_dosta.sh
#
# Harvest the dosta data from all of the OOI Coastal Endurance moorings. Data
# sets include telemetered and recovered host and instrument data. Data is
# downloaded from OOI Net and reworked to create a cleaner and more consistent
# set of files named and organized by the mooring, mooring sub-location, data
# delivery method and deployment.
#
# C. Wingard, 2020-08-21 -- Initial code

# set the base directory python command for all subsequent processing
. $(dirname $CONDA_EXE)/../etc/profile.d/conda.sh
conda activate ooi
PYTHON="python -m ooi_data_explorations.uncabled.process_dosta"

### CE01ISSM ###
BASE_FLAGS="-s CE01ISSM -n RID16 -sn 03-DOSTAD000"
BASE_FILE="${HOME}/ooidata/m2m/ce01issm/nsif/dosta/ce01issm.nsif.dosta"
for i in $(seq -f "%02g" 11 12); do
    $PYTHON $BASE_FLAGS -mt telemetered -st dosta_abcdjm_ctdbp_dcl_instrument -t ctdbp -dp $i -o "$BASE_FILE.deploy$i.telemetered.dosta_abcdjm_ctdbp_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st dosta_abcdjm_ctdbp_dcl_instrument_recovered -t ctdbp -dp $i -o "$BASE_FILE.deploy$i.recovered_host.dosta_abcdjm_ctdbp_dcl_instrument_recovered.nc"
    $PYTHON $BASE_FLAGS -mt recovered_inst -st dosta_abcdjm_ctdbp_instrument -t ctdbp -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.dosta_abcdjm_ctdbp_instrument.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st dosta_abcdjm_ctdbp_dcl_instrument -t ctdbp -dp 13 -o "$BASE_FILE.deploy13.telemetered.dosta_abcdjm_ctdbp_dcl_instrument.nc"

BASE_FLAGS="-s CE01ISSM -n MFD37 -sn 03-DOSTAD000"
BASE_FILE="${HOME}/ooidata/m2m/ce01issm/seafloor/dosta/ce01issm.seafloor.dosta"
for i in $(seq -f "%02g" 11 12); do
    $PYTHON $BASE_FLAGS -mt telemetered -st dosta_abcdjm_ctdbp_dcl_instrument -t ctdbp -dp $i -o "$BASE_FILE.deploy$i.telemetered.dosta_abcdjm_ctdbp_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st dosta_abcdjm_ctdbp_dcl_instrument_recovered -t ctdbp -dp $i -o "$BASE_FILE.deploy$i.recovered_host.dosta_abcdjm_ctdbp_dcl_instrument_recovered.nc"
    $PYTHON $BASE_FLAGS -mt recovered_inst -st dosta_abcdjm_ctdbp_instrument -t ctdbp -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.dosta_abcdjm_ctdbp_instrument.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st dosta_abcdjm_ctdbp_dcl_instrument -t ctdbp -dp 13 -o "$BASE_FILE.deploy13.telemetered.dosta_abcdjm_ctdbp_dcl_instrument.nc"

### CE02SHSM ###
BASE_FLAGS="-s CE02SHSM -n RID27 -sn 04-DOSTAD000"
BASE_FILE="${HOME}/ooidata/m2m/ce02shsm/nsif/dosta/ce02shsm.nsif.dosta"
for i in $(seq -f "%02g" 9 10); do
    $PYTHON $BASE_FLAGS -mt telemetered -st dosta_abcdjm_dcl_instrument -t solo -dp $i -ba -o "$BASE_FILE.deploy$i.telemetered.dosta_abcdjm_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st dosta_abcdjm_dcl_instrument_recovered -t solo -dp $i -ba -o "$BASE_FILE.deploy$i.recovered_host.dosta_abcdjm_dcl_instrument_recovered.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st dosta_abcdjm_dcl_instrument -t solo -dp 11 -ba -o "$BASE_FILE.deploy11.telemetered.dosta_abcdjm_dcl_instrument.nc"

### CE04OSSM ###
BASE_FLAGS="-s CE04OSSM -n RID27 -sn 04-DOSTAD000"
BASE_FILE="${HOME}/ooidata/m2m/ce04ossm/nsif/dosta/ce04ossm.nsif.dosta"
for i in $(seq -f "%02g" 8 9); do
    $PYTHON $BASE_FLAGS -mt telemetered -st dosta_abcdjm_dcl_instrument -t solo -dp $i -ba -o "$BASE_FILE.deploy$i.telemetered.dosta_abcdjm_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st dosta_abcdjm_dcl_instrument_recovered -t solo -dp $i -ba -o "$BASE_FILE.deploy$i.recovered_host.dosta_abcdjm_dcl_instrument_recovered.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st dosta_abcdjm_dcl_instrument -t solo -dp 10 -ba -o "$BASE_FILE.deploy10.telemetered.dosta_abcdjm_dcl_instrument.nc"

### CE06ISSM ###
BASE_FLAGS="-s CE06ISSM -n RID16 -sn 03-DOSTAD000"
BASE_FILE="${HOME}/ooidata/m2m/ce06issm/nsif/dosta/ce06issm.nsif.dosta"
for i in $(seq -f "%02g" 10 11); do
    $PYTHON $BASE_FLAGS -mt telemetered -st dosta_abcdjm_ctdbp_dcl_instrument -t ctdbp -dp $i -o "$BASE_FILE.deploy$i.telemetered.dosta_abcdjm_ctdbp_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st dosta_abcdjm_ctdbp_dcl_instrument_recovered -t ctdbp -dp $i -o "$BASE_FILE.deploy$i.recovered_host.dosta_abcdjm_ctdbp_dcl_instrument_recovered.nc"
    $PYTHON $BASE_FLAGS -mt recovered_inst -st dosta_abcdjm_ctdbp_instrument -t ctdbp -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.dosta_abcdjm_ctdbp_instrument.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st dosta_abcdjm_ctdbp_dcl_instrument -t ctdbp -dp 12 -o "$BASE_FILE.deploy12.telemetered.dosta_abcdjm_ctdbp_dcl_instrument.nc"

BASE_FLAGS="-s CE06ISSM -n MFD37 -sn 03-DOSTAD000"
BASE_FILE="${HOME}/ooidata/m2m/ce06issm/seafloor/dosta/ce06issm.seafloor.dosta"
for i in $(seq -f "%02g" 10 11); do
    $PYTHON $BASE_FLAGS -mt telemetered -st dosta_abcdjm_ctdbp_dcl_instrument -t ctdbp -dp $i -o "$BASE_FILE.deploy$i.telemetered.dosta_abcdjm_ctdbp_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st dosta_abcdjm_ctdbp_dcl_instrument_recovered -t ctdbp -dp $i -o "$BASE_FILE.deploy$i.recovered_host.dosta_abcdjm_ctdbp_dcl_instrument_recovered.nc"
    $PYTHON $BASE_FLAGS -mt recovered_inst -st dosta_abcdjm_ctdbp_instrument -t ctdbp -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.dosta_abcdjm_ctdbp_instrument.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st dosta_abcdjm_ctdbp_dcl_instrument -t ctdbp -dp 12 -o "$BASE_FILE.deploy12.telemetered.dosta_abcdjm_ctdbp_dcl_instrument.nc"

### CE07SHSM ###
BASE_FLAGS="-s CE07SHSM -n RID27 -sn 04-DOSTAD000"
BASE_FILE="${HOME}/ooidata/m2m/ce07shsm/nsif/dosta/ce07shsm.nsif.dosta"
for i in $(seq -f "%02g" 9 10); do
    $PYTHON $BASE_FLAGS -mt telemetered -st dosta_abcdjm_dcl_instrument -t solo -dp $i -ba -o "$BASE_FILE.deploy$i.telemetered.dosta_abcdjm_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st dosta_abcdjm_dcl_instrument_recovered -t solo -dp $i -ba -o "$BASE_FILE.deploy$i.recovered_host.dosta_abcdjm_dcl_instrument_recovered.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st dosta_abcdjm_dcl_instrument -t solo -dp 11 -ba -o "$BASE_FILE.deploy11.telemetered.dosta_abcdjm_dcl_instrument.nc"

BASE_FLAGS="-s CE07SHSM -n MFD37 -sn 03-DOSTAD000"
BASE_FILE="${HOME}/ooidata/m2m/ce07shsm/seafloor/dosta/ce07shsm.seafloor.dosta"
for i in $(seq -f "%02g" 9 10); do
    $PYTHON $BASE_FLAGS -mt telemetered -st dosta_abcdjm_ctdbp_dcl_instrument -t ctdbp -dp $i -o "$BASE_FILE.deploy$i.telemetered.dosta_abcdjm_ctdbp_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st dosta_abcdjm_ctdbp_dcl_instrument_recovered -t ctdbp -dp $i -o "$BASE_FILE.deploy$i.recovered_host.dosta_abcdjm_ctdbp_dcl_instrument_recovered.nc"
    $PYTHON $BASE_FLAGS -mt recovered_inst -st dosta_abcdjm_ctdbp_instrument -t ctdbp -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.dosta_abcdjm_ctdbp_instrument.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st dosta_abcdjm_ctdbp_dcl_instrument -t ctdbp -dp 11 -o "$BASE_FILE.deploy11.telemetered.dosta_abcdjm_ctdbp_dcl_instrument.nc"

### CE09OSSM ###
BASE_FLAGS="-s CE09OSSM -n RID27 -sn 04-DOSTAD000"
BASE_FILE="${HOME}/ooidata/m2m/ce09ossm/nsif/dosta/ce09ossm.nsif.dosta"
for i in $(seq -f "%02g" 9 10); do
    $PYTHON $BASE_FLAGS -mt telemetered -st dosta_abcdjm_dcl_instrument -t solo -dp $i -ba -o "$BASE_FILE.deploy$i.telemetered.dosta_abcdjm_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st dosta_abcdjm_dcl_instrument_recovered -t solo -dp $i -ba -o "$BASE_FILE.deploy$i.recovered_host.dosta_abcdjm_dcl_instrument_recovered.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st dosta_abcdjm_dcl_instrument -t solo -dp 11 -ba -o "$BASE_FILE.deploy11.telemetered.dosta_abcdjm_dcl_instrument.nc"

BASE_FLAGS="-s CE09OSSM -n MFD37 -sn 03-CTDBPE000"
BASE_FILE="${HOME}/ooidata/m2m/ce09ossm/seafloor/dosta/ce09ossm.seafloor.dosta"
for i in $(seq -f "%02g" 9 10); do
    $PYTHON $BASE_FLAGS -mt telemetered -st dosta_abcdjm_ctdbp_dcl_instrument -t ctdbp -dp $i -o "$BASE_FILE.deploy$i.telemetered.dosta_abcdjm_ctdbp_dcl_instrument.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st dosta_abcdjm_ctdbp_dcl_instrument_recovered -t ctdbp -dp $i -o "$BASE_FILE.deploy$i.recovered_host.dosta_abcdjm_ctdbp_dcl_instrument_recovered.nc"
    $PYTHON $BASE_FLAGS -mt recovered_inst -st dosta_abcdjm_ctdbp_instrument -t ctdbp -dp $i -o "$BASE_FILE.deploy$i.recovered_inst.dosta_abcdjm_ctdbp_instrument.nc"
done
# Current deployment
$PYTHON $BASE_FLAGS -mt telemetered -st dosta_abcdjm_ctdbp_dcl_instrument -t ctdbp -dp 11 -o "$BASE_FILE.deploy11.telemetered.dosta_abcdjm_ctdbp_dcl_instrument.nc"
