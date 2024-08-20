import os
from datetime import datetime

import numpy as np
from fdchp_utils import particles_to_pandas, read_file, read_file_to_pandas
from process_fdchp import process_fdchp

PA=0
EA=0
IS=1
if PA:
    lat=40.1334        #Pioneer NES
    Rwaves=[-0.75, 0, -5]
elif IS:
    lat=59.9337        #Irminger Sea
    Rwaves=[-0.75, 0, -6] 
else:
    lat=44.6393        #Endurance
    Rwaves=[-0.5, -0.5, -5]

U=[]
uw=[]
sigH=[]
data = []
file_path = '/c/Data/flux_calculations/FDCHPTester'
for incr in np.arange(24):
    #****************************************
    # Read in raw data
    #****************************************
    start = datetime.now()
    if incr <= 9:
        filename = os.path.join(file_path, 'fdchp_20200818_0'+ str(incr)+ '0200.dat')
    else:
        filename = os.path.join(file_path, 'fdchp_20200818_'+ str(incr)+ '0200.dat')

    print("Processing file: {}".format(filename))
    data = read_file(filename)

    raw_data = particles_to_pandas(data)
    # raw_data = read_file_to_pandas(filename)

    data_readin = datetime.now()
    
    #*****************************************
    # Compute flux data
    #*****************************************   

    fluxes, Uearth, waveheight = process_fdchp(raw_data, lat, Rwaves, file_suffix=str(incr))
    processed = datetime.now()
    print("File processing finished. Read-in time: {}, Process time: {}".format(data_readin-start, processed-data_readin))
    uw = uw + [-fluxes[0]]        # Fluxes: uw vw wT
    U = U + [Uearth]              # Wind speed relative to earth  
    sigH = sigH + [waveheight]    # Significant wave height

uw = np.array(uw)
U = np.array(U)
sigH = np.array(sigH)

print("Fin")