import numpy as np
from examples.fdchp.fdchp_utils import (
    alignwind,
    calculate_and_write_metrics,
    despikesimple,
    get_euler_angles,
    get_platform_vel_pos,
    gravity_from_latitude,
    heave_calc,
    sonic,
    write_metrics,
)
from scipy.signal import detrend, filtfilt


def process_fdchp(raw_data, latitude, anemometer_relative_position, tc1=20, tcwave=30, despike=True, despikecompass=False, flux_filepath=None):
    # JBE 06/29 JW 11Aug2014
    # 1.5 JW 3Sep2014
    # 2.0 JBE 06/12/2021
    # 3.0 JBE 05/11/2022
    # 3.5 JBE 01/02/2023
    # 4.0 JKP 2024-07-26 (pythonification)

    """
    tc1 is the cutoff period that determines where the pitch and roll Euler angles are 
    determined by the tilt of the accelerometers at low frequencies versus the 
    integrated angular rates at high frequencies.  This is known as complementary filtering

    tcwave is the same cutoff period as tc1, but for the wave measurements, which requires a large cutoff period
    """
    
    # Note: raw_data is expected to be in NWU coordinates, which is the default returned by the fdchp_utils.particles_to_pandas() function
    dt = 0.1 # Sampling period for FDCHP; 10 Hz
    fs = 1/dt # Sampling frequency for Windmaster

    # Define constants for filters
    tc2=tc1
    Rvec=np.zeros(3)             # Distance vector between the sonic anemometer and motion sensors
    Rvec[2]=0.85   # z offset
    FREEZING_POINT_K = 273.15

    G = gravity_from_latitude(latitude)

    version_number=4.0
    status_val = 1 #uint32

    #JBE Redefine files for 10 Hz and tc1=12 or 15 or 20.  Use digits(16) and
    #vpa(ahiwaves)
    if tc1==12:        #JBE Redefine filters for 10 Hz and tc1=12 with ludo=1
        ahi=[1.000000000000000,  -3.869797539975555,   5.617802044587569,  -3.625896801659086,   0.877898078061702]
        bhi=[0.9, 36962154017745,  -3.747848616070978,   5.621772924106467,  -3.747848616070978,   0.936962154017745]
    elif tc1==15:           #JBE Redefine filters for 10 Hz and tc1=15 with ludo=1
        ahi=[1.000000000000000,  -3.895833876325376,   5.692892648957240,  -3.698121672490409,   0.901065298297354]
        bhi=[0.949244593504399,  -3.796978374017595,   5.695467561026392,  -3.796978374017595,   0.949244593504399]
    else:
        ahi=[1.000000000000000,  -3.921872982100935,   5.768656400578301,  -3.771625137827138,   0.924842488052324]
        bhi=[0.961687313034919,  -3.846749252139674,   5.770123878209512,  -3.846749252139674,   0.961687313034919]

    #JBE Redefine files for 10 Hz and fcwaves=1/40;
    if tcwave==40:
        ahiwaves = [1.000000000000000,  -3.960935321365416,  5.883567180614652, -3.884319737084527,  0.961687926819144]
        bhiwaves = [0.980656885367734,  -3.922627541470935,  5.883941312206403, -3.922627541470935,  0.980656885367734]
    else:
        #JBE Redefine files for 10 Hz and fcwaves=1/30;
        ahiwaves = [1.0000000000000000, -3.947914166208924,  5.845094680183927, -3.846426389902994,  0.9492460297427443]
        bhiwaves = [0.9742925791274119, -3.897170316509647,  5.845755474764471, -3.897170316509647,  0.9742925791274119]

    num_datapoints=len(raw_data)
    mean_time = raw_data['time'].mean()

    #Some indices to be used later
    edge = int(np.fix( 20 * fs))
    tot=int(len(raw_data)-edge*2)
    start_index = int(edge)
    end_index = int(start_index+tot-1)
    
    if tot <= 0:
        # The datset is too short
        print("Unexpectedly short dataset in process_fdchp with shape: {}".format(raw_data.shape))
        return None
    
    #UNITS Velocities seem to be stored as cm/s, and must be converted to m/s
    #TODO: is this right for the case where this already represents a temperature?
    # This is the sonic temperature computed from speed of sound
    sos=raw_data['fdchp_speed_of_sound_sonic'].to_numpy()*0.01 
    if np.nanmedian(sos) < 50:
        Tv = sos
    else:
        # We convert SoS to temperature T1=C1^2/403 
        # print("adjusting speed of sound")
        Tv = (sos**2)/403.0 - FREEZING_POINT_K

    if despike:
        Tv = despikesimple(Tv)

    # sonic velocities
     # These are the 3-axis sonic wind velocities
     # Apparently stored as cm/s on disk; converted here to m/s.
    sonics=raw_data[['fdchp_wind_x', 'fdchp_wind_y', 'fdchp_wind_z']].to_numpy()*0.01

    if despike:                                # Remove obvious spikes in the data
        sonics = despikesimple(sonics)


    Tavg = np.nanmean(Tv)
    Tstd = np.nanstd(Tv)
    Tvmax = np.nanmax(Tv)
    Tvmin = np.nanmin(Tv)

    #***********************************************
    # Deal with the compass 
    # JBE Fill in bad points first units and signs 
    # UNITS-Roll, pitch and yaw are in radians
    #***********************************************

    compass=raw_data['fdchp_heading']  # heading
    roll=raw_data['fdchp_roll']
    pitch=raw_data['fdchp_pitch'] 


    gx=np.cos(compass)
    gy=np.sin(compass)
    compcos = np.mean(gx)
    compsin = np.mean(gy)
    compavg = np.arctan2(compsin,compcos)
    if compavg < 0:
        compavg=compavg + 2.0*np.pi

    if despikecompass:
        gx = despikesimple(gx)
        gy = despikesimple(gy)
        gsmooth=np.arctan2(gy,gx)
        gyro=gsmooth
    else:
        gyro=compass

    #***********************************************
    # Then the angular rates
    # UNITS - Rates are in radian/sec
    #***********************************************
    # This will hold the angular rates
    ang_rates = raw_data[['fdchp_x_ang_rate','fdchp_y_ang_rate','fdchp_z_ang_rate']].to_numpy()
    if despike:
        ang_rates = despikesimple(ang_rates)

    #***********************************************
    # Then the accelerations
    #***********************************************
    # This will hold the accelerations
    platform_accelerations = raw_data[['fdchp_x_accel_g', 'fdchp_y_accel_g', 'fdchp_z_accel_g']].to_numpy()
    if despike:
        platform_accelerations = despikesimple(platform_accelerations)

    #***************************************************************
    # We no longer force the mean accelerometers to equal gravity
    # UNITS#8 convert platform_accelerations accels to m/s^2
    #***************************************************************
    platform_accelerations=platform_accelerations*G

    #*********************************************
    #  Waves first
    #*********************************************
    # rpy = euler angles; dr = ang_rates rotated
    #TODO: change rpy and dr shapes to be 12000x3 instead of 3x12000
    rpy,dr = get_euler_angles(ahiwaves,bhiwaves,fs,platform_accelerations,ang_rates,gyro,G)  # euler angles are right-handed
    #[acc, uvwplatwave, xyzplat] = get_platform_vel_pos(bhiwaves,ahiwaves,fs,platform_accelerations,rpy,G);
    uvwplatwave,xyzplat = heave_calc(dr,rpy,platform_accelerations,fs,bhiwaves,ahiwaves,anemometer_relative_position,G)

    #*********************************************
    #  Then sonics
    #*********************************************
    euler,dr = get_euler_angles(ahi,bhi,fs,platform_accelerations,ang_rates,gyro,G)    # euler angles are right-handed
    acc, uvwplat, nope = get_platform_vel_pos(bhi,ahi,fs,platform_accelerations,euler,G)
    uvw,uvwr,uvwrot = sonic(sonics,dr,euler,uvwplat,Rvec)

    UVW = uvw[start_index:end_index, :]
    Ts = Tv[start_index:end_index]

    heave = detrend(xyzplat[start_index:end_index, 2]) #z-axis to find heave
    waveheight = 4*np.std(heave)

    U = np.mean(UVW[:, 0])
    V = np.mean(UVW[:, 1])
    W = np.mean(UVW[:, 2])

    # UNITS - wdir is in radians
    wdir = np.arctan2(V,-U)
    if wdir < 0:
        wdir = wdir + 2*np.pi

    u, alpha, beta = alignwind(UVW)

    # Relative to Earth
    Uearth=np.mean(u[0:tot, 0])

    for col in np.arange(u.shape[-1]):    # TODO: check axis
        u[:, col] = detrend(u[:, col])
    
    # Ts=detrend(Ts)
    
    Ts = detrend(Ts)

    fluxes = np.zeros(3)
    uwavg = np.mean(u[0:tot, 2]*u[0:tot, 0]) # Along-Wind component of momentum flux
    vwavg = np.mean(u[0:tot, 2]*u[0:tot, 1]) # Cross-Wind component of momentum flux
    wTavg = np.mean(u[0:tot, 2]*Ts) # Heat flux (vertical velocity times temp)
    #JBE Fixed a bug here by replacing mean with std

    Ucorstd = np.std(u[0:tot, 0])
    Vcorstd = np.mean(u[0:tot, 1]*u[0:tot, 1])
    Wcorstd = np.mean(u[0:tot, 2]*u[0:tot, 2])

    fluxes[0]=uwavg
    fluxes[1]=vwavg
    fluxes[2]=wTavg

    Uavg=np.mean(sonics[:, 0])
    Vavg=np.mean(sonics[:, 1])
    Wavg=np.mean(sonics[:, 2])
    rdir=np.arctan2(Vavg, Uavg) #TODO: not used?

    Umax=np.max(sonics[:, 0])
    Vmax=np.max(sonics[:, 1])
    Wmax=np.max(sonics[:, 2])
    Umin=np.min(sonics[:, 0])
    Vmin=np.min(sonics[:, 1])
    Wmin=np.min(sonics[:, 2])

    #JBE  This is not the standard deviation. 
    #Ustd = mean(abs(sonics(1,1:L)-Uavg));
    Ustd = np.std(sonics[:, 0])
    Vstd = np.std(sonics[:, 1])
    Wstd = np.std(sonics[:, 2])
    
    compstd=np.std(np.unwrap(compass))
    compmin=compass.min()
    compmax=compass.max()
    
    gx=np.cos(compass)
    gy=np.sin(compass)
    compcos = np.mean(gx)
    compsin = np.mean(gy)
    compavg = np.arctan2(compsin, compcos)
    
    mean_ang_rate=np.mean(ang_rates, axis=0)
    std_ang_rate=np.std(ang_rates, axis=0)

    Phiavg = mean_ang_rate[0]
    Thetaavg = mean_ang_rate[1]
    Psiavg = mean_ang_rate[2]
    Phistd = std_ang_rate[0]
    Thetastd = std_ang_rate[1]
    Psistd = std_ang_rate[2]
    Phimax = np.max(ang_rates[:, 0])
    Thetamax = np.max(ang_rates[:, 1])
    Psimax = np.max(ang_rates[:, 2])
    Phimin = np.min(ang_rates[:, 0])
    Thetamin = np.min(ang_rates[:, 1])
    Psimin = np.min(ang_rates[:, 2])
    
    gcomp = np.mean(platform_accelerations, axis=0)/G

    Axavg = gcomp[0]
    Ayavg = gcomp[1]
    Azavg = gcomp[2]
    Axmax=np.max(platform_accelerations[:, 0])/G
    Aymax=np.max(platform_accelerations[:, 1])/G
    Azmax=np.max(platform_accelerations[:, 2])/G
    Axmin=np.min(platform_accelerations[:, 0])/G
    Aymin=np.min(platform_accelerations[:, 1])/G
    Azmin=np.min(platform_accelerations[:, 2])/G

    gstd = np.std(platform_accelerations, axis=0)/G
    Axstd = gstd[0]/G
    Aystd = gstd[1]/G
    Azstd = gstd[2]/G
    
    # UNITS - Cal mean roll and pitch in radians
    Pitch = np.mean(pitch)
    Roll = np.mean(roll)
    Pitchstd = np.std(pitch)
    Rollstd = np.std(roll)
    Pitchmax = np.max(pitch)
    Rollmax = np.max(roll)
    Pitchmin = np.min(pitch)
    Rollmin = np.min(roll)

    version_number = 0
    status_val = 0

    flux_out = [
        mean_time,
        version_number,
        status_val,
        Uavg,
        Vavg,
        Wavg,
        Tavg,
        Ustd,
        Vstd,
        Wstd,
        Tstd,
        Umax,
        Vmax,
        Wmax,
        Tvmax,
        Umin,
        Vmin,
        Wmin,
        Tvmin,
        Axavg, #20
        Ayavg,
        Azavg,
        Axstd,
        Aystd,
        Azstd,
        Axmax,
        Aymax,
        Azmax,
        Axmin,
        Aymin,
        Azmin,
        Phiavg, #32
        Thetaavg,
        Psiavg,
        Phistd,
        Thetastd,
        Psistd,
        Phimax,
        Thetamax,
        Psimax,
        Phimin,
        Thetamin,
        Psimin,
        compavg,
        Pitch,
        Roll,
        compstd,
        Pitchstd,
        Rollstd,
        compmax,
        Pitchmax,
        Rollmax,
        compmin,
        Pitchmin,
        Rollmin,
        U,
        V,
        W,
        Ucorstd,
        Vcorstd,
        Wcorstd,
        Uearth,
        uwavg,
        vwavg,
        wTavg,
        waveheight
    ]

    file_path = flux_filepath if flux_filepath else "fluxes"
    write_metrics(flux_out, file_path=file_path)
    # calculate_and_write_metrics(mean_time, sonics, ang_rates, platform_accelerations, compass, roll, pitch, waveheight, fluxes)
    
    return (fluxes, Uearth, waveheight)



import numpy as np
from examples.fdchp.fdchp_utils import (
    alignwind,
    despikesimple,
    get_euler_angles,
    get_platform_vel_pos,
    gravity_from_latitude,
    heave_calc,
    sonic,
    write_metrics,
)
from scipy.signal import detrend


def process_fdchp_xarray(dataset, latitude, anemometer_relative_position, tc1=20, tcwave=30, despike=True, despikecompass=False, flux_filepath=None):
    # JBE 06/29 JW 11Aug2014
    # 1.5 JW 3Sep2014
    # 2.0 JBE 06/12/2021
    # 3.0 JBE 05/11/2022
    # 3.5 JBE 01/02/2023
    # 4.0 JKP 2024-07-26 (pythonification)
    # 5.0 JKP 2025-04-03 (xarray)

    """
    tc1 is the cutoff period that determines where the pitch and roll Euler angles are 
    determined by the tilt of the accelerometers at low frequencies versus the 
    integrated angular rates at high frequencies.  This is known as complementary filtering

    tcwave is the same cutoff period as tc1, but for the wave measurements, which requires a large cutoff period
    """
    
    # Note: dataset is expected to be in NWU coordinates, which is the default returned by the fdchp_utils.particles_to_pandas() function
    dt = 0.1 # Sampling period for FDCHP; 10 Hz
    fs = 1/dt # Sampling frequency for Windmaster

    # Define constants for filters
    Rvec=np.zeros(3)    # Distance vector between the sonic anemometer and motion sensors
    Rvec[2]=0.85        # z offset
    FREEZING_POINT_K = 273.15

    G = gravity_from_latitude(latitude)

    version_number=4.0
    status_val = 1 #uint32

    #JBE Redefine files for 10 Hz and tc1=12 or 15 or 20.  Use digits(16) and
    #vpa(ahiwaves)
    if tc1==12:        #JBE Redefine filters for 10 Hz and tc1=12 with ludo=1
        ahi=[1.000000000000000,  -3.869797539975555,   5.617802044587569,  -3.625896801659086,   0.877898078061702]
        bhi=[0.9, 36962154017745,  -3.747848616070978,   5.621772924106467,  -3.747848616070978,   0.936962154017745]
    elif tc1==15:           #JBE Redefine filters for 10 Hz and tc1=15 with ludo=1
        ahi=[1.000000000000000,  -3.895833876325376,   5.692892648957240,  -3.698121672490409,   0.901065298297354]
        bhi=[0.949244593504399,  -3.796978374017595,   5.695467561026392,  -3.796978374017595,   0.949244593504399]
    else:
        ahi=[1.000000000000000,  -3.921872982100935,   5.768656400578301,  -3.771625137827138,   0.924842488052324]
        bhi=[0.961687313034919,  -3.846749252139674,   5.770123878209512,  -3.846749252139674,   0.961687313034919]

    #JBE Redefine files for 10 Hz and fcwaves=1/40;
    if tcwave==40:
        ahiwaves = [1.000000000000000,  -3.960935321365416,  5.883567180614652, -3.884319737084527,  0.961687926819144]
        bhiwaves = [0.980656885367734,  -3.922627541470935,  5.883941312206403, -3.922627541470935,  0.980656885367734]
    else:
        #JBE Redefine files for 10 Hz and fcwaves=1/30;
        ahiwaves = [1.0000000000000000, -3.947914166208924,  5.845094680183927, -3.846426389902994,  0.9492460297427443]
        bhiwaves = [0.9742925791274119, -3.897170316509647,  5.845755474764471, -3.897170316509647,  0.9742925791274119]

    num_datapoints=len(dataset)
    mean_time = dataset['time'].mean()

    #Some indices to be used later
    edge = int(np.fix( 20 * fs))
    tot=int(len(dataset['time'])-edge*2)
    start_index = int(edge)
    end_index = int(start_index+tot-1)
    
    if tot <= 0:
        # The datset is too short
        print("Unexpectedly short dataset in process_fdchp with length: {}".format(dataset['time'].shape))
        return None
    
    #UNITS Velocities seem to be stored as cm/s, and must be converted to m/s
    #TODO: is this right for the case where this already represents a temperature?
    # This is the sonic temperature computed from speed of sound
    sos=dataset['fdchp_speed_of_sound_sonic'].to_numpy().T*0.01 
    if np.nanmedian(sos) < 50:
        Tv = sos
    else:
        # We convert SoS to temperature T1=C1^2/403 
        # print("adjusting speed of sound")
        Tv = (sos**2)/403.0 - FREEZING_POINT_K

    if despike:
        Tv = despikesimple(Tv)

    # sonic velocities
     # These are the 3-axis sonic wind velocities
     # Apparently stored as cm/s on disk; converted here to m/s.
    sonics=dataset[['fdchp_wind_x', 'fdchp_wind_y', 'fdchp_wind_z']].to_array().to_numpy().T*0.01

    if despike:                                # Remove obvious spikes in the data
        sonics = despikesimple(sonics)


    Tavg = np.nanmean(Tv)
    Tstd = np.nanstd(Tv)
    Tvmax = np.nanmax(Tv)
    Tvmin = np.nanmin(Tv)

    #***********************************************
    # Deal with the compass 
    # JBE Fill in bad points first units and signs 
    # UNITS-Roll, pitch and yaw are in radians
    #***********************************************

    compass=dataset['fdchp_heading'].T  # heading
    roll=dataset['fdchp_roll'].T
    pitch=dataset['fdchp_pitch'].T 


    gx=np.cos(compass)
    gy=np.sin(compass)
    compcos = np.mean(gx)
    compsin = np.mean(gy)
    compavg = np.arctan2(compsin,compcos)
    if compavg < 0:
        compavg=compavg + 2.0*np.pi

    if despikecompass:
        gx = despikesimple(gx)
        gy = despikesimple(gy)
        gsmooth=np.arctan2(gy,gx)
        gyro=gsmooth
    else:
        gyro=compass

    #***********************************************
    # Then the angular rates
    # UNITS - Rates are in radian/sec
    #***********************************************
    # This will hold the angular rates
    ang_rates = dataset[['fdchp_x_ang_rate','fdchp_y_ang_rate','fdchp_z_ang_rate']].to_array().to_numpy().T
    if despike:
        ang_rates = despikesimple(ang_rates)

    #***********************************************
    # Then the accelerations
    #***********************************************
    # This will hold the accelerations
    platform_accelerations = dataset[['fdchp_x_accel_g', 'fdchp_y_accel_g', 'fdchp_z_accel_g']].to_array().to_numpy().T
    if despike:
        platform_accelerations = despikesimple(platform_accelerations)

    #***************************************************************
    # We no longer force the mean accelerometers to equal gravity
    # UNITS#8 convert platform_accelerations accels to m/s^2
    #***************************************************************
    platform_accelerations=platform_accelerations*G

    #*********************************************
    #  Waves first
    #*********************************************
    # rpy = euler angles; dr = ang_rates rotated
    #TODO: change rpy and dr shapes to be 12000x3 instead of 3x12000
    rpy,dr = get_euler_angles(ahiwaves,bhiwaves,fs,platform_accelerations,ang_rates,gyro,G)  # euler angles are right-handed
    #[acc, uvwplatwave, xyzplat] = get_platform_vel_pos(bhiwaves,ahiwaves,fs,platform_accelerations,rpy,G);
    uvwplatwave,xyzplat = heave_calc(dr,rpy,platform_accelerations,fs,bhiwaves,ahiwaves,anemometer_relative_position,G)

    #*********************************************
    #  Then sonics
    #*********************************************
    euler,dr = get_euler_angles(ahi,bhi,fs,platform_accelerations,ang_rates,gyro,G)    # euler angles are right-handed
    acc, uvwplat, nope = get_platform_vel_pos(bhi,ahi,fs,platform_accelerations,euler,G)
    uvw,uvwr,uvwrot = sonic(sonics,dr,euler,uvwplat,Rvec)

    UVW = uvw[start_index:end_index, :]
    Ts = Tv[start_index:end_index]

    heave = detrend(xyzplat[start_index:end_index, 2]) #z-axis to find heave
    waveheight = 4*np.std(heave)

    U = np.mean(UVW[:, 0])
    V = np.mean(UVW[:, 1])
    W = np.mean(UVW[:, 2])

    # UNITS - wdir is in radians
    wdir = np.arctan2(V,-U)
    if wdir < 0:
        wdir = wdir + 2*np.pi

    u, alpha, beta = alignwind(UVW)

    # Relative to Earth
    Uearth=np.mean(u[0:tot, 0])

    for col in np.arange(u.shape[-1]):    # TODO: check axis
        u[:, col] = detrend(u[:, col])
    
    # Ts=detrend(Ts)
    
    Ts = detrend(Ts)

    fluxes = np.zeros(3)
    uwavg = np.mean(u[0:tot, 2]*u[0:tot, 0]) # Along-Wind component of momentum flux
    vwavg = np.mean(u[0:tot, 2]*u[0:tot, 1]) # Cross-Wind component of momentum flux
    wTavg = np.mean(u[0:tot, 2]*Ts) # Heat flux (vertical velocity times temp)

    fluxes[0]=uwavg
    fluxes[1]=vwavg
    fluxes[2]=wTavg

    return (fluxes, Uearth, waveheight)
