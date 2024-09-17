import os
from datetime import datetime
from glob import glob

import numpy as np
from fdchp_utils import particles_to_pandas, read_file, read_file_to_pandas
from process_fdchp import process_fdchp
from plot_fdchp import plot_x_velocity_vs_wind_speed, plot_wave_height_vs_wind_speed

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
# file_path = '/c/Data/flux_calculations/FDCHPTester'

# directory = '/home/jovyan/ooi/uncabled/GI01SUMO/R00008/instruments/dcl12/FDCHP_sn144904/D202105'
directory = '/home/jovyan/ooi/uncabled/GI01SUMO/R00008/instruments/dcl12/FDCHP_sn144904/D202110'

files = glob(os.path.join(directory, '*.dat'))
output_filepath = "fluxes"
output_filename = "fluxes{}"

if not os.path.exists(output_filepath):
    os.makedirs(output_filepath)
print("Processing {} files.".format(len(files)))
incr=0
print("Start time: {}".format(datetime.now()))
# for incr in np.arange(24):
for filename in files:
    #****************************************
    # Read in raw data
    #****************************************
    start = datetime.now()
    # if incr <= 9:
    #     filename = os.path.join(file_path, 'fdchp_20200818_0'+ str(incr)+ '0200.dat')
    # else:
    #     filename = os.path.join(file_path, 'fdchp_20200818_'+ str(incr)+ '0200.dat')

    print("Processing file {}: {}".format(incr + 1, filename))
    try:
        data = read_file(filename)
    except Exception as e:
        print("Error reading in file: {}".format(e))
        continue
              
    raw_data = particles_to_pandas(data)

    # raw_data = read_file_to_pandas(filename)

    data_readin = datetime.now()
    
    #*****************************************
    # Compute flux data
    #*****************************************   
    try:
        fluxes, Uearth, waveheight = process_fdchp(raw_data, lat, Rwaves, flux_filepath=os.path.join(output_filepath, output_filename.format(incr)))
    except Exception as e:
        # Error processing data; probably too few datapoints
        print("Error processing fdchp dataset: {}".format(e))
        continue
    processed = datetime.now()
    print("File processing finished. Read-in time: {}, Process time: {}".format(data_readin-start, processed-data_readin))
    uw = uw + [-fluxes[0]]        # Fluxes: uw vw wT
    U = U + [Uearth]              # Wind speed relative to earth  
    sigH = sigH + [waveheight]    # Significant wave height
    incr += 1

uw = np.array(uw)
U = np.array(U)
sigH = np.array(sigH)


plot_x_velocity_vs_wind_speed(U, uw)

plot_wave_height_vs_wind_speed(U, sigH)

print("End time: {}".format(datetime.now()))
print("Fin")