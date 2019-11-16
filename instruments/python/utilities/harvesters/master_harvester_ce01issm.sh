#!/usr/bin/env bash

cd /cygdrive/c/Users/cwingard/Documents/GitHub/ooi-data-explorations
SITE="CE01ISSM"
DEPLOY=11
METHOD="telemetered"

### BUOY
LOCATE="buoy"
NODE="SBD17"

# ctdbp
INSTRMT="ctdbp"
SENSOR="06-CTDBPC000"
STREAM="ctdbp_cdef_dcl_instrument"
FILE="${SITE,,}/$LOCATE/$INSTRMT/${SITE,,}.$LOCATE.$INSTRMT.deploy$DEPLOY.$METHOD.$STREAM.nc"
FLAGS="-s $SITE -n $NODE -sn $SENSOR -mt $METHOD -st $STREAM -dp $DEPLOY"
python -m instruments.python.$INSTRMT.request_$INSTRMT $FLAGS -o $FILE
# FLORT
INSTRMT="flort"
SENSOR="06-FLORTD000"
STREAM="flort_sample"
FILE="${SITE,,}/$LOCATE/$INSTRMT/${SITE,,}.$LOCATE.$INSTRMT.deploy$DEPLOY.$METHOD.$STREAM.nc"
FLAGS="-s $SITE -n $NODE -sn $SENSOR -mt $METHOD -st $STREAM -dp $DEPLOY"
python -m instruments.python.$INSTRMT.request_$INSTRMT $FLAGS -o $FILE
# VELPT
# MOPAK

### NSIF
# ctdbp
# DOSTA
# FLORT
# NUTNR
# OPTAA
# PC02W
# PHSEN
# SPKIR
# VELPT

### MFN
# ADCPT
# ctdbp
# DOSTA
# OPTAA
# PC02W
# PHSEN
# PRESF
# VEL3D
# ZPLSC
