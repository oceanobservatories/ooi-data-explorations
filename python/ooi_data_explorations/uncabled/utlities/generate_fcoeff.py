#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import numpy as np
import os
import sys


def create_fcoeff(dspec_file, fcoeff_file, site_depth):
    """
    Creates the non-directional spectra parameters needed to recreate the RDI
    WavesMon v3 FCoeff file(s) from the ADCP wave directional spectra (DSpec)
    file(s). The FCoeff file is what the OOI CI system is expecting for the
    ADCPT-M WAVSS non-directional spectra data products. With the upgrade to
    the WavesMon v4 processing utility, these files are no longer generated.
    This utility will recreate them so the data can be added to the OOI data
    store. Note, this does not recreate 1:1 the original FCoeff files, as the
    original files were created using a different method than the one used
    here. The results are very similar, however, and as the original algorithm
    is not available, this is the best that can be done using a different
    method described in the CDIP document ADCPwaves2CDIP.pdf referenced below.

    :param dspec_file: Directional spectra file (DSpec) to be processed
    :param fcoeff_file: FCoeff file to be created from the DSpec file
    :param site_depth: Site depth used to set an upper cutoff frequency for the
        non-directional spectra parameters
    :return: non-directional spectra and mean wave direction

    reference: https://cordc.ucsd.edu/about/docs/ADCPwaves2CDIP.pdf
    """
    # calculate the upper cutoff frequency (see references above)
    hm = site_depth / 2.0  # half the site depth, mid-water depth
    theta = np.deg2rad(20.0)  # 20 degrees, the angle of the ADCP beams converted to radians
    d = np.sqrt(2.0) * hm * np.sin(theta)  # distance between adjacent transducers at mid-water depth

    k = 2.0 * np.pi / (2 * d)  # wave number
    w2 = 9.81 * k * np.tanh(k * site_depth)  # wave frequency squared
    upper_cutoff = np.sqrt(w2) / (2.0 * np.pi)  # upper cutoff frequency

    # now load the DSpec data file
    with open(dspec_file, 'r') as f:
        data = f.readlines()

    # parse the header
    header = data[2].split()
    ndir = int(header[1])  # number of directions
    nfreq = int(header[4])  # number of frequencies
    cutoff = float(header[-1])  # cutoff frequency
    header = data[4].replace('(', ' ').replace(')', '').split()
    band_width = float(header[4])  # frequency bandwidth
    band_start = float(header[-1])  # starting frequency, mid-band

    # The WavesMon software reports a cutoff frequency, which we always want to use unless it is greater than the
    # upper cutoff frequency calculated above. If it is, then we want to use the upper cutoff frequency instead.
    if cutoff > upper_cutoff:
        cutoff = upper_cutoff

    # calculate the frequency array and a mask for the cutoff frequencies
    start = band_start + (band_width / (nfreq / 2))
    frequency = np.arange(start, start + (band_width * nfreq), band_width)
    mask = (frequency < 0.01) | (frequency > cutoff)

    # check to see if the DSpec data has been screened by the ADCP software. If it has, then there is no data to
    # process, and we can return None for the non-directional spectra and mean wave direction.
    if data[6] == 'Screened\n':
        return None, None, False

    # if not, parse the frequency by direction data into an array
    freq = [[np.int64(x) for x in f.split()] for f in data[6:] if f.split()]
    freq = np.array(freq) * 1e-6  # convert from mm^2/Hz to m^2/Hz

    # calculate the non-directional spectra and convert from mm^2/Hz to m^2/Hz
    non_dir = np.sum(freq, axis=1) / ndir

    # calculate the mean wave direction
    directions = np.deg2rad(np.arange(0, 360, 360 / ndir))
    mean_dir = np.arctan2(np.sum(freq * np.sin(directions), axis=1), np.sum(freq * np.cos(directions), axis=1))

    # calculate the fourier coefficients (a1, a2, b1 and b2)
    a0 = np.sum(freq, axis=1)
    n = a0 > 0.0
    a1 = np.zeros(a0.shape)
    a1[n] = 1 / a0[n] * np.sum(freq[n] * np.cos(directions), axis=1)
    a2 = np.zeros(a0.shape)
    a2[n] = 1 / a0[n] * np.sum(freq[n] * np.cos(2 * directions), axis=1)
    b1 = np.zeros(a0.shape)
    b1[n] = 1 / a0[n] * np.sum(freq[n] * np.sin(directions), axis=1)
    b2 = np.zeros(a0.shape)
    b2[n] = 1 / a0[n] * np.sum(freq[n] * np.sin(2 * directions), axis=1)

    # can check the results above by calculating the mean direction using the fourier coefficients. the mean direction
    # calculated this way should be equal to the mean direction calculated above to at least the 6th decimal place.
    check = np.isclose(mean_dir, np.arctan2(b1, a1), rtol=1e-06, atol=1e-08)

    # now reset the mean directions to degrees and make sure they are between 0 and 360
    mean_dir = np.rad2deg(mean_dir)
    n = mean_dir < 0.0
    mean_dir[n] = mean_dir[n] + 360.0

    # apply the frequency cutoff mask to the non-directional spectra, mean wave direction and fourier coefficients
    non_dir[mask] = 0.0
    mean_dir[mask] = 0.0
    a1[mask] = 0.0
    a2[mask] = 0.0
    b1[mask] = 0.0
    b2[mask] = 0.0

    # create the FCoeffYYMMDDhhmm.txt file
    with open(fcoeff_file, 'w+') as f:
        # write the header
        f.write('\n')
        f.write(' % Fourier Coefficients\n')
        f.write('% 9 Fields and {:d} Frequencies\n'.format(nfreq))
        f.write('% Frequency(Hz), Band width(Hz), Energy density(m^2/Hz), Direction (deg), '
                'A1, B1, A2, B2, Check Factor\n')
        f.write('% Frequency Bands are {:.6f} Hz wide (first frequency band is centered '
                'at {:.8f})\n'.format(band_width, band_start))
        # now write the data
        for i in range(nfreq):
            f.write(' {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} -9999.9\n'.format(frequency[i],
                                                                                                band_width,
                                                                                                non_dir[i], mean_dir[i],
                                                                                                a1[i].round(6) + 0.0,
                                                                                                b1[i].round(6) + 0.0,
                                                                                                a2[i].round(6) + 0.0,
                                                                                                b2[i].round(6) + 0.0))

    # return the non-directional spectra and mean wave direction arrays for cross-checking purposes
    return non_dir, mean_dir, check


def log9_to_fcoeff(log9_file):
    """
    Use the data in the Log 9 formatted bulk wave statistics file to create the
    non-directional spectra parameters needed to recreate the RDI WavesMon v3
    FCoeff file(s) from the ADCP wave directional spectra (DSpec) file(s).

    :param log9_file: Log 9 formatted bulk wave statistics file
    :return: None
    """
    # load the log9 bulk wave statistics file
    base_dir = os.path.dirname(log9_file)
    with open(log9_file, 'r') as f:
        data = f.readlines()

    # for each burst measurement, use the date and time to create a file name for the FCoeff file, pull out the
    # site depth from the water level and the peak mean direction and then create the FCoeff file.
    for line in data:
        # split the data and extract the date and time information, water level and peak mean direction
        burst = line.split(',')
        year = burst[1]
        month = burst[2]
        day = burst[3]
        hour = burst[4]
        minute = burst[5]
        water_level = int(burst[17]) / 1000.0 - 1.0  # convert from mm to m and subtract 1 m for the ADCP offset
        peak_mean_direction = int(burst[10])

        # create the file names for the DSpec and FCoeff files and then create the FCoeff file
        dspec_file = 'DSpec' + year + month + day + hour + minute + '.txt'
        dspec_file = os.path.join(base_dir, dspec_file)
        fcoeff_file = 'FCoeff' + year + month + day + hour + minute + '.txt'
        fcoeff_file = os.path.join(base_dir, fcoeff_file)
        non_dir, mean_dir, check = create_fcoeff(dspec_file, fcoeff_file, water_level)

        # check to make sure the two different methods for calculating the mean direction are equivalent
        if non_dir is not None:
            #m = non_dir.argmax()
            #pd = np.abs(mean_dir[m] - peak_mean_direction) / peak_mean_direction
            #if pd > 0.10:
            #    print('Mean directions do not align for {:s}. The peak mean direction ({:d}) and the '
            #          'non-directional estimate ({:.1f}) differ by more than 10%. Marking as failed.'
            #          .format(fcoeff_file, peak_mean_direction, mean_dir[m]))
            #    os.rename(fcoeff_file, fcoeff_file + '.fail')
            if not check.all():
                print('Internal consistency check failure for {:s}'.format(fcoeff_file))
                os.rename(fcoeff_file, fcoeff_file + '.fail')
        else:
            print('DSpec file {:s} was screened by WaveMon processing, skipping.'.format(dspec_file))


def inputs(argv=None):
    """
    Sets the main input arguments that are needed to run the script.

    :param argv: command line arguments
    :return: input arguments
    """
    if argv is None:
        argv = sys.argv[1:]

    # assign input arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--log9', dest='log9_file', type=str, required=True,
                        help='Log 9 formatted bulk wave statistics file')

    # parse the input arguments and create a parser object
    args = parser.parse_args(argv)

    return args


def main(argv=None):
    """
    Main entry point for the script. Parses the command line arguments and
    calls the appropriate function to create the FCoeff file(s). The FCoeff
    file(s) are what the OOI CI system is expecting for the ADCPT-M WAVSS
    non-directional spectra data products.

    :param argv: command line arguments
    :return: None
    """
    # parse the command line arguments
    args = inputs(argv)

    # call the function to create the FCoeff file(s)
    log9_to_fcoeff(args.log9_file)


if __name__ == '__main__':
    main()
