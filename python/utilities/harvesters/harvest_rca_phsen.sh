#!/usr/bin/env bash
#
# harvest_rca_phsen.sh
#
# Harvest the phsen data from all of the OOI Regional Cabled Array systems. Data 
# is downloaded from OOI Net and reworked to create a cleaner and more consistent
# set of files named and organized by the platform, platform sub-location, data
# delivery method and deployment.
#
# C. Wingard, 2019-12-12 -- Initial code

# set the base directory python command for all subsequent processing
cd /home/ooiuser/code/ooi-data-explorations
PYTHON="python -m python.ooi_data_explorations.cabled.request_phsen"

### CE02SHBP ###
BASE_FLAGS="-s CE02SHBP -n LJ01D -sn 10-PHSEND103"
BASE_FILE="ce02shbp/seafloor/phsen/ce02shbp.seafloor.phsen"
for i in $(seq -f "%02g" 1 6); do
	$PYTHON $BASE_FLAGS -mt streamed -st phsen_data_record -dp $i -o "$BASE_FILE.deploy$i.streamed.phsen_data_record.nc"
done

### CE04OSBP ###
BASE_FLAGS="-s CE04OSBP -n LJ01C -sn 10-PHSEND107"
BASE_FILE="ce04osbp/seafloor/phsen/ce04osbp.seafloor.phsen"
for i in $(seq -f "%02g" 1 6); do
  $PYTHON $BASE_FLAGS -mt streamed -st phsen_data_record -dp $i -o "$BASE_FILE.deploy$i.streamed.phsen_data_record.nc"
done

### CE04OSPS ###
BASE_FLAGS="-s CE04OSPS -n PC01B -sn 4B-PHSENA106"
BASE_FILE="ce04osps/200m/phsen/ce04osps.200m.phsen"
for i in $(seq -f "%02g" 1 6); do
  $PYTHON $BASE_FLAGS -mt streamed -st phsen_data_record -dp $i -o "$BASE_FILE.deploy$i.streamed.phsen_data_record.nc"
done

BASE_FLAGS="-s CE04OSPS -n SF01B -sn 2B-PHSENA108"
BASE_FILE="ce04osps/profiler/phsen/ce04osps.profiler.phsen"
for i in $(seq -f "%02g" 1 6); do
  $PYTHON $BASE_FLAGS -mt streamed -st phsen_data_record -dp $i -o "$BASE_FILE.deploy$i.streamed.phsen_data_record.nc"
done

### RS01SBPS ###
BASE_FLAGS="-s RS01SBPS -n PC01A -sn 4B-PHSENA102"
BASE_FILE="rs01sbps/200m/phsen/rs01sbps.200m.phsen"
for i in $(seq -f "%02g" 1 6); do
  $PYTHON $BASE_FLAGS -mt streamed -st phsen_data_record -dp $i -o "$BASE_FILE.deploy$i.streamed.phsen_data_record.nc"
done

BASE_FLAGS="-s RS01SBPS -n SF01A -sn 2D-PHSENA101"
BASE_FILE="rs01sbps/profiler/phsen/rs01sbps.profiler.phsen"
for i in $(seq -f "%02g" 1 7); do
  $PYTHON $BASE_FLAGS -mt streamed -st phsen_data_record -dp $i -o "$BASE_FILE.deploy$i.streamed.phsen_data_record.nc"
done

### RS03AXPS ###
BASE_FLAGS="-s RS03AXPS -n PC03A -sn 4B-PHSENA302"
BASE_FILE="rs03axps/200m/phsen/rs03axps.200m.phsen"
for i in $(seq -f "%02g" 1 5); do
  $PYTHON $BASE_FLAGS -mt streamed -st phsen_data_record -dp $i -o "$BASE_FILE.deploy$i.streamed.phsen_data_record.nc"
done

BASE_FLAGS="-s RS03AXPS -n SF03A -sn 2D-PHSENA301"
BASE_FILE="rs03axps/profiler/phsen/rs03axps.profiler.phsen"
for i in $(seq -f "%02g" 1 5); do
  $PYTHON $BASE_FLAGS -mt streamed -st phsen_data_record -dp $i -o "$BASE_FILE.deploy$i.streamed.phsen_data_record.nc"
done
