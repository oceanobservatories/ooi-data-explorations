import datetime

import numpy as np
import pandas as pd
from mi.dataset.parser.fdchp_a import FdchpADataParticle, FdchpAParser
from scipy import integrate, interpolate
from scipy.signal import detrend, filtfilt

PADLENGTH = 12 # (3*(np.max([len(bhi), len(ahi)]) - 1) == 12)

def exception_handler(exception):
    print("Exception!".format(exception))
    
def read_file(file_path):
    data = []
    with open(file_path, 'rb') as input:
        parser = FdchpAParser(input, exception_handler)
        # parser.parse_file()
        particle = parser.get_records()
        
        while particle:
            data.append(particle[0])
            particle = parser.get_records()
    return data

def read_file_to_pandas(file_path):
    data = []
    with open(file_path, 'rb') as input:
        parser = FdchpAParser(input, exception_handler)
        # parser.parse_file()
        particle = parser.get_records()
        
        while particle:
            data.append({ value['value_id']: [value['value']] for value in particle[0].generate_dict()['values']})
            particle = parser.get_records()
    if data:
        df = pd.DataFrame(data)
    return df

def particles_to_pandas(particles, convert_to_nwu=True):
    df = None
    for particle in particles:
        values = particle.generate_dict()['values']
        data_dict = { value['value_id']: [value['value']] for value in values}
        entry = pd.DataFrame.from_dict(data_dict)
        df = pd.concat([df,entry])
        
    if df is not None:
        df['time'] = df.apply(lambda row: datetime.datetime( int(row.year), int(row.month), int(row.day), int(row.hour), int(row.minute), int(row.second), int(row.millisecond)*1000), axis=1)
        
        if convert_to_nwu:
            # Convert IMU from North East Down coordinate system to North West Up coordinate system to match Sonic.
            # Note that we can leave the x-axis variable as measured.

            df['fdchp_heading'] = -df['fdchp_heading']   # z heading(yaw) counter clockwise
            df['fdchp_pitch'] = -df['fdchp_pitch']   # y pitch east to west
            df['fdchp_y_ang_rate'] = -df['fdchp_y_ang_rate'] # angular rate around y-axis
            df['fdchp_z_ang_rate'] = -df['fdchp_z_ang_rate'] # angular rate around z-axis 
            df['fdchp_y_accel_g'] = -df['fdchp_y_accel_g'] # linear acceleration along y-axis
            df['fdchp_z_accel_g'] = -df['fdchp_z_accel_g'] # linear acceleration along z-axis (positve up)
        
    return df
    

def despikesimple(Y, exclusion_stdevs = 4, iterations = 3):
    # Remove Outliers by excluding points outside median + exclusion_stdevs* stdev
    # from an interpolation, then interpolating over the full dataset again
    # Y is numpy.array
    # iterations is number of passes
    # exclusion_stdevs is the number of standard deviations beyond which to exclude data

    # TODO: fix this function - seems not to work like the MATLAB version
    data_shape = Y.shape
    rows = data_shape[0]

    if len(data_shape) == 1:
        # one-dimensional
        cols = 1
    else:
        cols = data_shape[1]
    t = np.arange(rows)

    def handle_interp(rows, data):
        median = np.nanmedian(data)
        stdev = np.nanstd(data)
        good_data = (data < median + exclusion_stdevs*stdev) & (data > median - exclusion_stdevs*stdev)
        acceptable_points = data[good_data]
        k = np.sum(good_data)
        # if k != len(data):
        #     print("despikesimple despiking {} points".format(len(data) - k))
        if k == len(rows):
            return data# don't bother iterating when there aren't any points to exclude
        if k > 0:
            interp = interpolate.interp1d(rows[good_data], acceptable_points, fill_value='extrapolate')
            data[:] = interp(rows)

        return data

    for iter in range(iterations):
        if cols == 1:
            X = Y[:]
        else:
            for col in range(cols):
                X = Y[:, col]
                Y[:, col] = handle_interp(t, X)

    return Y
                

def gravity_from_latitude(latitude):
    # Get gravitational acceleration g based on latitude
    if latitude < 0 or latitude > 90:
        latitude = 45.0

    gamma = 9.7803267715
    c1 = 0.0052790414
    c2 = 0.0000232718
    c3 = 0.0000001262
    c4 = 0.0000000007

    phi = latitude*np.pi/180.0
    x = np.sin(phi)
    g = gamma*(1+c1*x**2+c2*x**4+c3*x**6+c4*x**8)

    return g


def rotate_buoy_to_world(input_angles, euler_angles, iflag = False):
    # Rotate from buoy frame to world frame

    sinp  = np.sin(euler_angles[:, 0])
    cosp  = np.cos(euler_angles[:, 0])
    sint  = np.sin(euler_angles[:, 1])
    cost  = np.cos(euler_angles[:, 1])
    sinps = np.sin(euler_angles[:, 2])
    cosps = np.cos(euler_angles[:, 2])

    up = input_angles[:, 0]
    vp = input_angles[:, 1]
    wp = input_angles[:, 2]

    if iflag: #  from xyz to x'y'z'
        u = up*cost*cosps                   + vp*cost*sinps                   - wp*sint
        v = up*(sinp*sint*cosps-cosp*sinps) + vp*(sinp*sint*sinps+cosp*cosps) + wp*(cost*sinp)
        w = up*(cosp*sint*cosps+sinp*sinps) + vp*(cosp*sint*sinps-sinp*cosps) + wp*(cost*cosp)

    else: # from x'y'z' to xyz
        u = up*cost*cosps + vp*(sinp*sint*cosps-cosp*sinps) + wp*(cosp*sint*cosps+sinp*sinps)
        v = up*cost*sinps + vp*(sinp*sint*sinps+cosp*cosps) + wp*(cosp*sint*sinps-sinp*cosps)
        w = up*(-sint)    + vp*(cost*sinp)                  + wp*(cost*cosp)

    return np.array([u, v, w]).T

def sonic(sonics, omegam, euler, uvwplat, R):
    """
    Function from EDDYCORR toolbox
    CORRECT SONIC ANEMOMETER COMPONENTS FOR PLATFORM MOTION AND ORIENTATION.
    INPUTS:
       sonics     - row of integers corre to sonic numbers which are to be 
          corrected
       omegam     - (3XN) measured angular rate 'vector' in platform frame
       euler      - (3XN) array of euler angles (phi, theta, psi)
       uvwplat    - (3XN) array of platform velocities (output from accels_.m)
       R          - Distance vector between the sonic anemometer and motion sensors

    OUTPUTS:
       uvw        - (MXN) array of corrected sonic anemometer components, in the 
                        fixed earth reference frame  (North-West-up)
    CALCULATE WINDS IN EARTH BASED FRAME. THE ANGULAR VELOCITY IS CALCULATED AS
    THE CROSS PRODUCT BETWEEN THE ANGULAR RATE VECTOR AND POSITION VECTOR.THE 
    MEASURED AND ANGULAR VELOCITIES ARE IN THE PLATFORM FRAME AND MUST BE
    ROTATED INTO THE EARTH FRAME. THE PLATFORM VELOCITY IS ALREADY IN THE EARTH
    FRAME (FROM ACCELS.M), SO HERE THEY CAN JUST BE ADDED.
     UVW =  MEASURED VELOCITY + ANGULAR RATE INDUCED VELOCITIES + 
      INTEGRATED ACCELEROMETERS  
  """

    Rvec = np.array(R * np.ones(omegam.shape))
    uvwrot = np.cross(omegam,Rvec)

    uvw  = rotate_buoy_to_world(sonics + uvwrot, euler) + uvwplat
    uvwr = rotate_buoy_to_world(sonics + uvwrot, euler)
    return (uvw, uvwr, uvwrot)



def get_angle_update_matrix(angular_rates, angles):
    """
    Function from EDDYCORR toolbox
    This function computes the angular update matrix
    as described in Edson et al. (1998) and Thwaites
    (1995) page 50.
    """
    p  = angles[:, 0]
    t  = angles[:, 1]
    ps = angles[:, 2]

    up = angular_rates[:, 0]
    vp = angular_rates[:, 1]
    wp = angular_rates[:, 2]

    u = up  + vp*np.sin(p)*np.tan(t) + wp*np.cos(p)*np.tan(t)
    v =  0  + vp*np.cos(p)           - wp*np.sin(p)
    w =  0  + vp*np.sin(p)/np.cos(t) + wp*np.cos(p)/np.cos(t)

    return np.array([u, v, w]).T

def get_euler_angles(ahi, bhi, sampling_frequency, accelerations, angular_rates, gyro, gravity, iterations = 5):
    """
     Function from EDDYCORR toolbox

     Sept 2018  Made compass right-handed by changing:
      psi_slow= -gyro - filtfilt(bhi,ahi,-gyro) to
      psi_slow=  gyro - filtfilt(bhi,ahi, gyro)

     Sept 2000     Replaced integrations with cumtrapz function

     May 16 1997 - modified to remove the first estimate of the euler
        angles in the nonlinear euler angle update matrix, F^-1 matrix
        is approximated by the identity matrix. still uses trapezoidal
        intetgration

     INPUT:
        ahi,bhi - 5 element array of filter coefficients for filtfilt
        sampling_frequency    - sampling frequency
        accelerations  - (3xN) array of recalibrated linear accelerations: acc_x,acc_y,acc_z
        angular_rates - (3xN) array of recalibrated angular rates: rate_x, rate_y, rate_z
        gyro  - (1xN) array of gyro signal
        iterations   - number of interations

     OUTPUT:
        euler    - (3XN) array of the euler angles (phi, theta, psi) in radians.
        dr       - ???

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     THE ANGLES ARE  ESTIMATED FROM

     angle = slow_angle (from accelerometers) + fast_angle (integrated rate sensors)
    """
    # Unwrap compass
    gyro = np.unwrap(gyro)

    ang_rates = np.copy(angular_rates)
    # REMOVE MEAN FROM RATE SENSORS
    for col in np.arange(ang_rates.shape[-1]):    # TODO: check axis
        ang_rates[:, col] = detrend(ang_rates[:,col])

    # LOW FREQUENCY ANGLES FROM ACCELEROMETERS AND GYRO
    # SLOW ROLL FROM GRAVITY EFFECTS ON HORIZONTAL ACCELERATIONS. LOW PASS
    # FILTER SINCE HIGH FREQUENCY HORIZONTAL ACCELERATIONS MAY BE 'REAL'

    # PITCH
    normalized_x_acceleration = -accelerations[:,0]/gravity
    theta = np.minimum(normalized_x_acceleration, 1) # Cap accelerations by g; small angles
    theta = np.maximum(theta,-1)                    
    sensible_x_acc =np.nonzero(np.abs(normalized_x_acceleration) < 1)  # Get indices where x-acceleration is less than g
    theta[sensible_x_acc] = np.arcsin(normalized_x_acceleration[sensible_x_acc])
    #TODO: check out filtfilt some more. It looks like scipy default filtering may be different than MATLAB's:
    # https://dsp.stackexchange.com/questions/11466/differences-between-python-and-matlab-filtfilt-function
    theta_slow = theta - filtfilt(bhi, ahi, theta, padlen = PADLENGTH) # 3*(np.max([len(bhi), len(ahi)]) - 1)) # padlen set to match MATLAB default

    # ROLL
    normalized_y_acceleration = accelerations[:, 1]/gravity
    phi = np.minimum(normalized_y_acceleration, 1) # Cap accelerations by g; small angles
    phi = np.maximum(phi,-1)
    acc_y_grav_cos_theta_slow = normalized_y_acceleration/np.cos(theta_slow)
    sensible_y_acc = np.nonzero(np.abs(acc_y_grav_cos_theta_slow) < 1)
    phi[sensible_y_acc] = np.arcsin(acc_y_grav_cos_theta_slow[sensible_y_acc]) # Get indices where y-acceleration is less than g
    phi_slow = phi - filtfilt(bhi, ahi, phi, padlen = PADLENGTH) # 3*(np.max([len(bhi), len(ahi)]) - 1)) # padlen set to match MATLAB default

    # YAW
    # HERE, WE ESTIMATE THE SLOW HEADING. THE 'FAST HEADING' IS NOT NEEDED
    # FOR THE EULER ANGLE UPDATE MATRIX.
    # Using the microstrain compass for heading
    # THE COMPASS IS PASSED RIGHT HANDED SO IT CAN BE TREATED LIKE PITCH AND ROLL
    psi_slow = gyro - filtfilt(bhi, ahi, gyro, padlen = PADLENGTH) #3*(max([len(bhi), len(ahi)]) - 1)) # padlen set to match MATLAB default

    # USE SLOW ANGLES AS FIRST GUESS
    euler  = np.array([phi_slow, theta_slow, psi_slow]).T
    rates = get_angle_update_matrix(ang_rates, euler)
    # INTEGRATE AND FILTER ANGLE RATES, AND ADD TO SLOW ANGLES
    for i in np.arange(iterations):
        #TODO: check column/row ordering of phi,theta, psi from rates
        phi_int = integrate.cumulative_trapezoid(rates[:, 0], dx = 1/sampling_frequency, initial=0)
        phi = phi_slow + filtfilt(bhi, ahi, phi_int, padlen = PADLENGTH) # 3*(np.max([len(bhi), len(ahi)]) - 1)) # padlen set to match MATLAB default 
        theta_int = integrate.cumulative_trapezoid(rates[:, 1], dx = 1/sampling_frequency, initial=0)
        theta = theta_slow + filtfilt(bhi, ahi, theta_int, padlen = PADLENGTH) # 3*(np.max([len(bhi), len(ahi)]) - 1)) # padlen set to match MATLAB default
        psi_int = integrate.cumulative_trapezoid(rates[:, 2], dx = 1/sampling_frequency, initial=0)
        psi = psi_slow + filtfilt(bhi, ahi, psi_int, padlen = PADLENGTH) #3*(np.max([len(bhi), len(ahi)]) - 1)) # padlen set to match MATLAB default 

        euler  = np.array([phi, theta, psi]).T
        rates = get_angle_update_matrix(ang_rates, euler)
        for col in np.arange(rates.shape[-1]):
            rates[:, col] = detrend(rates[:, col])
        #TODO: could I replace the detrend above with this simple mean subtraction?                       
        # rates = rates - np.mean(rates, axis=0)

    return (euler, rates)

def heave_calc(omegam, euler, accm, sampling_frequency, bhi, ahi, R, gravity):
    """
    Function from EDDYCORR toolbox
    CORRECT SONIC ANEMOMETER COMPONENTS FOR PLATFORM MOTION AND ORIENTATION.
    INPUTS:
       omegam     - (3XN) measured angular rate 'vector' in platform frame
       euler      - (3XN) array of euler angles (phi, theta, psi)
       accm       - (3XN) array of platform accelerations
       R          - vector distance from motionPak to wave sensor

    OUTPUTS:
      platform_velocity - platform velocity at sensor location: [u, v, w]
      platform_disp - platform displacement at sensor location: [x, y, z]
    """

    Rvec = np.array(R * np.ones(accm.shape))
    uvwrot = np.cross(omegam, Rvec)
    uvwrot  = rotate_buoy_to_world(uvwrot, euler)

    acc = rotate_buoy_to_world(accm, euler) # first rotate
    acc[:, 2] = acc[:, 2] - gravity      # remove gravity from z-component of acceleration

    platform_velocity = np.zeros(acc.shape)
    platform_disp = np.zeros(acc.shape)
    
    for i in range(3):
        platform_velocity[:, i] = integrate.cumulative_trapezoid(acc[:, i], dx = 1/sampling_frequency, initial=0) + uvwrot[:, i]
        platform_velocity[:, i] = filtfilt(bhi, ahi, platform_velocity[:, i], padlen = PADLENGTH) # 3*(np.max([len(bhi), len(ahi)]) - 1)) # padlen set to match MATLAB default 
        # INTEGRATE AGAIN TO GET DISPLACEMENTS
        platform_disp[:, i] = integrate.cumulative_trapezoid(platform_velocity[:, i], dx = 1/sampling_frequency, initial=0)
        platform_disp[:, i] = filtfilt(bhi, ahi, platform_disp[:, i], padlen = PADLENGTH) #3*(np.max([len(bhi), len(ahi)]) - 1)) # padlen set to match MATLAB default 

    return (platform_velocity, platform_disp)

def get_platform_vel_pos(bhi, ahi, sampling_frequency, accm, euler, gravity):
    """
    Function from EDDYCORR toolbox
    2008          Allows filter to have different cutoff from angular filter
    Sept 2000     Replaced integrations with cumtrapz function
    Mar 3 1998	Redesigned the high pass filter (see below).
    Revised: June 10, 1997 - high pass filter with higher cutoff 
               frequency than previous version
    Calculates platform velocity and positions from accelerometer data.
    Integrates linear accelerations to get platform velocity
    and displacement. After each integration, signals are 
    high pass filtered to remove low frequency effects.
    INPUT
       bhigh,ahigh - high pass filter coefficients
       sampling_frequency      - sampling frequency
       accm    - calibrated linear accelerations (output from recal.m)
       euler   - (3xN) Euler angles phi, theta, psi 
    OUTPUT:
       acc     - (3XN) linear accelerations in FLIP/Earth reference 
       lin_velocities - (3XN) linear velocities at the point of motion measurement: [u, v, w]
       platform_disp - (3XN) platform displacements from mean position: [x, y, z]
    """
    acc = rotate_buoy_to_world(accm, euler) # first rotate
    acc[:, 2] = acc[:, 2] - gravity # remove gravity

    # INTEGRATE ACCELERATIONS TO GET PLATFORM VELOCITIES
    lin_velocities = np.zeros(acc.shape)
    platform_disp = np.zeros(acc.shape)
    for i in range(3):
        lin_velocities[:, i] = integrate.cumulative_trapezoid(acc[:, i], dx = 1/sampling_frequency, initial=0)
        lin_velocities[:, i] = filtfilt(bhi, ahi, lin_velocities[:, i], padlen = PADLENGTH)
        # INTEGRATE AGAIN TO GET DISPLACEMENTS
        platform_disp[:, i] = integrate.cumulative_trapezoid(lin_velocities[:, i], dx = 1/sampling_frequency, initial=0)
        platform_disp[:, i] = filtfilt(bhi, ahi, platform_disp[:, i], padlen = PADLENGTH)

    return (acc, lin_velocities, platform_disp)

def alignwind(U):
    # u,v,w are in platform coordinates;
    # rotate motion corrected velocities into mean wind
    # fluxes are relative to mean wind (for whole 20-minute segment)

    Ub = np.mean(U[:, 0])
    Vb = np.mean(U[:, 1])
    Wb = np.mean(U[:, 2])
    Sb = np.sqrt(Ub**2 + Vb**2)
    beta  = np.arctan2(Wb,Sb)
    alpha = np.arctan2(Vb,Ub)
    Ur =  U[:, 0]*np.cos(alpha)*np.cos(beta) + U[:, 1]*np.sin(alpha)*np.cos(beta) + U[:,2]*np.sin(beta)
    Vr = -U[:, 0]*np.sin(alpha) + U[:, 1]*np.cos(alpha)
    Wr = -U[:, 0]*np.cos(alpha)*np.sin(beta) - U[:, 1]*np.sin(alpha)*np.sin(beta) + U[:, 2]*np.cos(beta)

    u=np.zeros(U.shape)

    u[:, 0] = Ur
    u[:, 1] = Vr
    u[:, 2] = Wr

    beta  = beta*180/np.pi
    alpha = alpha*180/np.pi
    return [u, alpha, beta]

def calculate_and_write_metrics(mean_time, sonics, ang_rates, platform_accelerations, compass, roll, pitch, fluxes, waveheight, file_path='fluxes'):
    
    Uavg=np.mean(sonics[:, 0])
    Vavg=np.mean(sonics[:, 1])
    Wavg=np.mean(sonics[:, 2])
    rdir=np.arctan2(Vavg,Uavg) #TODO: not used?

    Umax=np.max(sonics[0,:])
    Vmax=np.max(sonics[1,:])
    Wmax=np.max(sonics[2,:])
    Umin=np.min(sonics[0,:])
    Vmin=np.min(sonics[1,:])
    Wmin=np.min(sonics[2,:])

    #JBE  This is not the standard deviation. 
    #Ustd = mean(abs(sonics(1,1:L)-Uavg));
    Ustd = np.std(sonics[0,:])
    Vstd = np.std(sonics[1,:])
    Wstd = np.std(sonics[2,:])
    
    compstd=np.std(np.unwrap(compass))
    compmin=compass.min()
    compmax=compass.max()
    
    gx=np.cos(compass)
    gy=np.sin(compass)
    compcos = np.mean(gx)
    compsin = np.mean(gy)
    compavg = np.arctan2(compsin,compcos)
    
    mean_ang_rate=np.mean(ang_rates, axis=0) #TODO: check axis
    std_ang_rate=np.std(ang_rates, axis=0) #TODO: check axis

    Phiavg = mean_ang_rate[0]
    Thetaavg = mean_ang_rate[1]
    Psiavg = mean_ang_rate[2]
    Phistd = std_ang_rate[0]
    Thetastd = std_ang_rate[1]
    Psistd = std_ang_rate[2]
    Phimax = np.max(ang_rates[0,:])
    Thetamax = np.max(ang_rates[1,:])
    Psimax = np.max(ang_rates[2,:])
    Phimin = np.min(ang_rates[0,:])
    Thetamin = np.min(ang_rates[1,:])
    Psimin = np.min(ang_rates[2,:])
    
    gcomp = np.mean(platform_accelerations, axis=0)

    Axavg = gcomp[0]
    Ayavg = gcomp[1]
    Azavg = gcomp[2]
    Axmax=np.max(platform_accelerations[:,0])
    Aymax=np.max(platform_accelerations[:, 1])
    Azmax=np.max(platform_accelerations[:, 2])
    Axmin=np.min(platform_accelerations[:, 0])
    Aymin=np.min(platform_accelerations[:, 1])
    Azmin=np.min(platform_accelerations[:, 2])

    gstd = np.std(platform_accelerations, axis=0)
    Axstd = gstd[0]
    Aystd = gstd[1]
    Azstd = gstd[2]
    
    # UNITS - Cal mean roll and pitch in radians
    Pitch = np.mean(pitch)
    Roll = np.mean(roll)
    Pitchstd = np.std(pitch)
    Rollstd = np.std(roll)
    Pitchmax = np.max(pitch)
    Rollmax = np.max(roll)
    Pitchmin = np.min(pitch)
    Rollmin = np.min(roll)
                 
    with open(file_path, 'w') as outfile:
        outfile.write(str(fluxes))
        
def write_metrics(metrics_list, file_path='fluxes', rounding=2):
    with open(file_path, 'w') as outfile:
        for metric in metrics_list:
            met = np.round(metric, rounding) if np.isscalar(metric) else metric
            outfile.write(str(met) + ", ")
