#!/usr/bin/env bash
#
# harvest_ce_vel3d.sh
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
PYTHON="python -m ooi_data_explorations.uncabled.process_vel3d"

### CE01ISSM ###
BASE_FLAGS="-s CE01ISSM -n MFD35 -sn 01-VEL3DD000"
BASE_FILE="${HOME}/ooidata/m2m/ce01issm/mfn/vel3d/ce01issm.mfn.vel3d"
for i in $(seq -f "%02g" 1 19); do
    $PYTHON $BASE_FLAGS -mt telemetered -st vel3d_cd_dcl_velocity_data -ba -db $i -o "$BASE_FILE.deploy$i.telemetered.vel3d_cd_dcl_velocity_data.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st vel3d_cd_dcl_velocity_data_recovered -ba -db $i -o "$BASE_FILE.deploy$i.recovered_host.vel3d_cd_dcl_velocity_data_recovered.nc"
done

### CE06ISSM ###
BASE_FLAGS="-s CE06ISSM -n MFD35 -sn 01-VEL3DD000"
BASE_FILE="${HOME}/ooidata/m2m/ce06issm/mfn/vel3d/ce06issm.mfn.vel3d"
for i in $(seq -f "%02g" 1 18); do
    $PYTHON $BASE_FLAGS -mt telemetered -st vel3d_cd_dcl_velocity_data -ba -db $i -o "$BASE_FILE.deploy$i.telemetered.vel3d_cd_dcl_velocity_data.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st vel3d_cd_dcl_velocity_data_recovered -ba -db $i -o "$BASE_FILE.deploy$i.recovered_host.vel3d_cd_dcl_velocity_data_recovered.nc"
done

### CE07SHSM ###
BASE_FLAGS="-s CE07SHSM -n MFD35 -sn 01-VEL3DD000"
BASE_FILE="${HOME}/ooidata/m2m/ce07shsm/mfn/vel3d/ce07shsm.mfn.vel3d"
for i in $(seq -f "%02g" 1 17); do
    $PYTHON $BASE_FLAGS -mt telemetered -st vel3d_cd_dcl_velocity_data -ba -db $i -o "$BASE_FILE.deploy$i.telemetered.vel3d_cd_dcl_velocity_data.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st vel3d_cd_dcl_velocity_data_recovered -ba -db $i -o "$BASE_FILE.deploy$i.recovered_host.vel3d_cd_dcl_velocity_data_recovered.nc"
done

### CE09OSSM ###
BASE_FLAGS="-s CE09OSSM -n MFD35 -sn 01-VEL3DD000"
BASE_FILE="${HOME}/ooidata/m2m/ce09ossm/mfn/vel3d/ce09ossm.mfn.vel3d"
for i in $(seq -f "%02g" 1 17); do
    $PYTHON $BASE_FLAGS -mt telemetered -st vel3d_cd_dcl_velocity_data -ba -db $i -o "$BASE_FILE.deploy$i.telemetered.vel3d_cd_dcl_velocity_data.nc"
    $PYTHON $BASE_FLAGS -mt recovered_host -st vel3d_cd_dcl_velocity_data_recovered -ba -db $i -o "$BASE_FILE.deploy$i.recovered_host.vel3d_cd_dcl_velocity_data_recovered.nc"
done
