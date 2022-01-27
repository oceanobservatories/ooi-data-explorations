#!/usr/bin/env bash
#
# harvest_ce_profiler_flort.sh
#
# Harvest the FLORT data from all of the OOI Coastal Endurance profiler
# moorings. Data sets include telemetered and recovered profiler data. Data is
# downloaded from OOINet and reworked to create a cleaner and more consistent
# set of files named and organized by the mooring, mooring sub-location, data
# delivery method and deployment.
#
# C. Wingard, 2019-07-22 -- Initial code

# set the base directory python command for all subsequent processing
. $(dirname $CONDA_EXE)/../etc/profile.d/conda.sh
conda activate ooi
PYTHON="python -m ooi_data_explorations.uncabled.process_flort"

### CE01ISSP ###
BASE_FLAGS="-s CE01ISSP -n SP001 -sn 08-FLORTJ000"
BASE_FILE="${HOME}/ooidata/m2m/ce01issp/flort/ce01issp.flort"
for i in $(seq -f "%02g" 1 17); do
    #$PYTHON $BASE_FLAGS -mt telemetered -st flort_sample -dp $i -o "$BASE_FILE.deploy$i.telemetered.flort_sample.nc"
    $PYTHON $BASE_FLAGS -mt recovered_cspp -st flort_sample -dp $i -o "$BASE_FILE.deploy$i.recovered_cspp.flort_sample.nc"
done

### CE02SHSP ###
BASE_FLAGS="-s CE02SHSP -n SP001 -sn 07-FLORTJ000"
BASE_FILE="${HOME}/ooidata/m2m/ce02shsp/flort/ce02shsp.flort"
for i in $(seq -f "%02g" 1 20); do
    #$PYTHON $BASE_FLAGS -mt telemetered -st flort_sample -dp $i -o "$BASE_FILE.deploy$i.telemetered.flort_sample.nc"
    $PYTHON $BASE_FLAGS -mt recovered_cspp -st flort_sample -dp $i -o "$BASE_FILE.deploy$i.recovered_cspp.flort_sample.nc"
done

### CE06ISSP ###
BASE_FLAGS="-s CE06ISSP -n SP001 -sn 08-FLORTJ000"
BASE_FILE="${HOME}/ooidata/m2m/ce06issp/flort/ce06issp.flort"
for i in $(seq -f "%02g" 1 13); do
    #$PYTHON $BASE_FLAGS -mt telemetered -st flort_sample -dp $i -o "$BASE_FILE.deploy$i.telemetered.flort_sample.nc"
    $PYTHON $BASE_FLAGS -mt recovered_cspp -st flort_sample -dp $i -o "$BASE_FILE.deploy$i.recovered_cspp.flort_sample.nc"
done

### CE07SHSP ###
BASE_FLAGS="-s CE07SHSP -n SP001 -sn 07-FLORTJ000"
BASE_FILE="${HOME}/ooidata/m2m/ce07shsp/flort/ce07shsp.flort"
for i in $(seq -f "%02g" 1 10); do
    #$PYTHON $BASE_FLAGS -mt telemetered -st flort_sample -dp $i -o "$BASE_FILE.deploy$i.telemetered.flort_sample.nc"
    $PYTHON $BASE_FLAGS -mt recovered_cspp -st flort_sample -dp $i -o "$BASE_FILE.deploy$i.recovered_cspp.flort_sample.nc"
done

### CE09OSPM ###
BASE_FLAGS="-s CE09OSPM -n WFP01 -sn 04-FLORTK000"
BASE_FILE="${HOME}/ooidata/m2m/ce09ospm/flort/ce09ospm.flort"
for i in $(seq -f "%02g" 1 15); do
    $PYTHON $BASE_FLAGS -mt telemetered -st flort_sample -dp $i -o "$BASE_FILE.deploy$i.telemetered.flort_sample.nc"
    $PYTHON $BASE_FLAGS -mt recovered_wfp -st flort_sample -dp $i -o "$BASE_FILE.deploy$i.recovered_wfp.flort_sample.nc"
done
