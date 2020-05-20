#!/usr/bin/env bash

cd /home/ooiuser/code/ooi-data-explorations
SITE="CE01ISSM"
DEPLOY=12
METHOD="telemetered"

### BUOY
LOCATE="buoy"
NODE="SBD17"

# CTDBP
INSTRMT="ctdbp"
SENSOR="06-CTDBPC000"
STREAM="ctdbp_cdef_dcl_instrument"
FILE="${SITE,,}/$LOCATE/$INSTRMT/${SITE,,}.$LOCATE.$INSTRMT.deploy$DEPLOY.$METHOD.$STREAM.nc"
FLAGS="-s $SITE -n $NODE -sn $SENSOR -mt $METHOD -st $STREAM -dp $DEPLOY"
python -m ooi_data_explorations.uncabled.request_$INSTRMT $FLAGS -o $FILE
# FLORT
INSTRMT="flort"
SENSOR="06-FLORTD000"
STREAM="flort_sample"
FILE="${SITE,,}/$LOCATE/$INSTRMT/${SITE,,}.$LOCATE.$INSTRMT.deploy$DEPLOY.$METHOD.$STREAM.nc"
FLAGS="-s $SITE -n $NODE -sn $SENSOR -mt $METHOD -st $STREAM -dp $DEPLOY"
python -m ooi_data_explorations.uncabled.request_$INSTRMT $FLAGS -o $FILE
# VELPT
# MOPAK

### NSIF
LOCATE="nsif"
NODE="RID16"

# CTDBP
INSTRMT="ctdbp"
SENSOR="03-CTDBPC000"
STREAM="ctdbp_cdef_dcl_instrument"
FILE="${SITE,,}/$LOCATE/$INSTRMT/${SITE,,}.$LOCATE.$INSTRMT.deploy$DEPLOY.$METHOD.$STREAM.nc"
FLAGS="-s $SITE -n $NODE -sn $SENSOR -mt $METHOD -st $STREAM -dp $DEPLOY"
python -m ooi_data_explorations.uncabled.request_$INSTRMT $FLAGS -o $FILE
# FLORT
INSTRMT="flort"
SENSOR="02-FLORTD000"
STREAM="flort_sample"
FILE="${SITE,,}/$LOCATE/$INSTRMT/${SITE,,}.$LOCATE.$INSTRMT.deploy$DEPLOY.$METHOD.$STREAM.nc"
FLAGS="-s $SITE -n $NODE -sn $SENSOR -mt $METHOD -st $STREAM -dp $DEPLOY"
python -m ooi_data_explorations.uncabled.request_$INSTRMT $FLAGS -o $FILE
# DOSTA
# FLORT
# NUTNR
# OPTAA
# PC02W
# PHSEN
INSTRMT="phsen"
SENSOR="06-PHSEND000"
STREAM="phsen_abcdef_dcl_instrument"
FILE="${SITE,,}/$LOCATE/$INSTRMT/${SITE,,}.$LOCATE.$INSTRMT.deploy$DEPLOY.$METHOD.$STREAM.nc"
FLAGS="-s $SITE -n $NODE -sn $SENSOR -mt $METHOD -st $STREAM -dp $DEPLOY"
python -m ooi_data_explorations.uncabled.request_$INSTRMT $FLAGS -o $FILE
# SPKIR
# VELPT

### MFN
LOCATE="seafloor"
NODE="MFD35"

# ADCPT
# PC02W
# PHSEN
INSTRMT="phsen"
SENSOR="06-PHSEND000"
STREAM="phsen_abcdef_dcl_instrument"
FILE="${SITE,,}/$LOCATE/$INSTRMT/${SITE,,}.$LOCATE.$INSTRMT.deploy$DEPLOY.$METHOD.$STREAM.nc"
FLAGS="-s $SITE -n $NODE -sn $SENSOR -mt $METHOD -st $STREAM -dp $DEPLOY"
python -m ooi_data_explorations.uncabled.request_$INSTRMT $FLAGS -o $FILE
# PRESF
# VEL3D

NODE="MFD37"
# CTDBP
INSTRMT="ctdbp"
SENSOR="03-CTDBPC000"
STREAM="ctdbp_cdef_dcl_instrument"
FILE="${SITE,,}/$LOCATE/$INSTRMT/${SITE,,}.$LOCATE.$INSTRMT.deploy$DEPLOY.$METHOD.$STREAM.nc"
FLAGS="-s $SITE -n $NODE -sn $SENSOR -mt $METHOD -st $STREAM -dp $DEPLOY"
python -m ooi_data_explorations.uncabled.request_$INSTRMT $FLAGS -o $FILE
# DOSTA
# OPTAA
# ZPLSC
