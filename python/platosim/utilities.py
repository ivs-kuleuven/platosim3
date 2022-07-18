#!/usr/bin/env python3

"""
This python module contains all general utilities that are commonly used
by the different codes within the PlatoSim and the PLATOnium repository.
"""

import os
import sys
import h5py
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import astropy.units as u
from astropy.coordinates import SkyCoord
from pylab import MaxNLocator
from colorama import Fore, Style
from scipy.ndimage import median_filter
from numba import njit

from platosim.simulation import Simulation
from platosim.referenceFrames import (skyToFocalPlaneCoordinates,
                                      gnomonicRadialDistanceFromOpticalAxis,
                                      getCCDandPixelCoordinates)


#==============================================================#
#                           FUNCTIONS                          #
#==============================================================#


def errorcode(API, message):
    """
    This function allows to colour code error messages within a code.
    """
    if API == 'software':
        print(Style.BRIGHT + Fore.BLUE + message + Style.RESET_ALL)
    if API == 'module':
        print(Style.BRIGHT + Fore.GREEN + message + Style.RESET_ALL)
    if API == 'message':
        print(Style.BRIGHT + message + Style.RESET_ALL)
    if API == 'warning':
        print(Style.BRIGHT + Fore.YELLOW + '[Warning]: ' + message + Style.RESET_ALL)
    if API == 'error':
        print(Style.BRIGHT + Fore.RED + '[Error]: ' + message + Style.RESET_ALL)
        exit()





def compilation(i, i_max, text=''):
    """
    This function print out a compilation-time-bar in the terminal.

    PARAMETERS
    ----------
    i : int
        Index of the current running job-loop
    i_max : int
        Index of the last running job-loop
    text : str
        Optional text written next to compilation-bar

    RETURN
    ------
        No parameters only a nice compilation bar to bash
    """

    # Running percentage calculation

    percent = (i + 1) / (i_max * 1.0) * 100

    # We here divide by 2 as the length of the bar is only 50 characters:

    bar = "[" + "-" * int(percent/2) + '>' + " " * (50-int(percent/2))+"] {}% {}"\
        .format(int(percent), text)

    # Print and clean-up

    sys.stdout.write(u"\u001b[1000D" + bar)
    sys.stdout.flush()





def normalize(signal, factor=1e6, length=-1):
    """
    This function normalize a signal with the option to factorize it.

    PARAMETERS
    ----------
    signal : narray
        Signal array that needs to be turned into relative signal.
    factor : int, float
        Factor to scale signal. Default is in ppm.
    length : int, float
        Option to normalize signal using only a smaller initial signal segment.
        Default is to use the entire signal sequence.

    RETURN
    ------
    relative_signal : narray
        Normalized relative signal is returned. Default unit in [ppm].
    """

    relative_signal = (signal / np.mean(signal[:int(length)]) - 1) * factor

    return relative_signal






@njit
def filter(signal, filt='median', carbox=144):
    """
    This utility makes the proper filter solution to a signal dataset.
    Notice: the carbox size here is twice what is default by numpy.

    PARAMETERS
    ----------
    filtType : string
        Filter is either: median or mean
    signalIn : narray
       Signal needed for processing
    numBox : int
        Integer used as car-box size of 1 hour: 3600s/25s = 144.

    RETURN
    ------
    signalOut : narray
        Filtered signal array
    """

    # Constants

    n     = carbox
    S     = signal.copy()     # Avoid overwritting the input signal
    S_new = np.zeros(len(S))  # Prepare forloop
    nzero = np.zeros(2*n+1)   # Optimization constant

    for i in range(len(S)-2*n):

        # Interval: d[n, 1+n, ... , N-1, N-n]

        if filt == 'median': S_new[n+i] = np.median(S[np.arange((n+i)-n, (n+i)+n+1)])
        if filt == 'mean':   S_new[n+i] = np.mean(S[np.arange((n+i)-n, (n+i)+n+1)])
        if filt == 'std':    S_new[n+i] = np.std(S[np.arange((n+i)-n, (n+i)+n+1)])

    for i in range(n):

        # Interval: d[-n, -(n-1), ... , n-1, n] - Low end of data

        low = nzero
        low[np.arange(n-i)] = S[0]*np.ones(n-i)
        low[-(n+1+i):] = S[np.arange(0, n+1+i)]

        if filt == 'median': S_new[i] = np.median(low)
        if filt == 'mean':   S_new[i] = np.mean(low)
        if filt == 'std':    S_new[i] = np.std(low)

        # Interval: d[N-n, N-(n-1), ... , N+(n-1), N+n] - High end of data

        high = nzero
        high[np.arange(n+1+i)] = S[np.arange(len(S)-(n+i+1), len(S))]
        high[-(n-i):]      = S[-1]*np.ones(n-i)

        if filt == 'median': S_new[len(S)-1-i] = np.median(high)
        if filt == 'mean':   S_new[len(S)-1-i] = np.mean(high)
        if filt == 'std':    S_new[len(S)-1-i] = np.std(high)

    return S_new






def passbandConversionV2P(V, Teff):
    """
    Coversion from Johnson-Cousin V magnitude to the PLATO passband.
    This filtersion is from Marchiori et al. (2019), Eq. 5 and 6, and is
    extracted using synthetic stellar spectra from the POLLOX database.
    NOTE valid for Teff = 4000-15000 K (hence not for M-dwarfs).

    PARAMETERS
    ----------
    V : float, narray
        Johnson-Cousin magnitude of star(s).
    Teff : float narray
        Effective temperature of star(s).

    RETURN
    ------
    P : float, narray
        The PLATO passband magnitude of star(s).
    """

    # The actual filtersion equation

    c = [1.184e-12, 4.526e-8, 5.805e-4, 2.449]     # Machiori et al. (2019)
    #c = [2.366e-12, 8.126e-08, -0.0009279, 3.499] # Fabio Fialho et al. in prep
    P  = c[0]*Teff**3 - c[1]*Teff**2 + c[2]*Teff - c[3] + V

    return P







def NSRphotonNoiseLimit(P, Ncam=24., Ntra=1., tdur=3600., camType='N'):
    """
    NSR estimate in the photon noise limit of bright stars. The stellar flux are
    calculated from the PLATO passband found by Marchiori et al. (2019).
    NOTE only valid for very bright stars (P < 11).

    PARAMETERS
    ----------
    P : float, narray
        The PLATO passband magnitude.
    Ncam : float, narray
        Number of telescope visibility. N-Cams (6, 12. 18, 24) or F-Cams (2).
    Ntra : float, narray
        Number of transits that can be co-added by phase-folding.
    tdur : float, narray
        Time duration over which the NSR is estimated. E.g., 3600s for 1h precision.
    camType : str
        Either the normal (N) or fast (F) cameras. Default is normal.

    RETURN
    ------
    NSR : float, narray
        NSR only valid for the photon noise limit.
    """

    # Choose cycle and exposure time for either the normal (N) or fast (F) cameras

    if camType == 'N':
        texp = 21.  # [s]
        tcyc = 25.  # [s]
    elif camType == 'F':
        texp = 2.1
        tcyc = 2.5

    # The P passband zero-point

    zp = 20.62

    # Flux of stars [e-/s]

    f = 10**(-0.4*(P-zp))

    # Observed total flux per exposure in ADU counts

    g = 900000/65535.  # [e-/ADU] Gain
    F = f * texp / g   # [ADU]

    # SNR from pure photon noise and NSR from uncorrelated noise.
    # Gaussian statistic gives sigma --> sigma/sqrt(N)

    SNR = np.sqrt(F * Ncam * Ntra * tdur/tcyc)
    NSR = 1/SNR * 1e6

    return NSR






def pdAddColumn(df, newCol, name):
    """
    Function to add a column to an exisiting pandas data frame.
    """
    df[name] = newCol
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    return df[cols]





def convertQuarterRange(dQ):
    """
    Small function that takes a string of numbers (here quarters)
    and split it up into readable float values used as real number
    ranges. If a single number is given, an quarter integer is 
    returned.
    """
    quarters = []
    for part in dQ.split(','):
        if '-' in part:
            # If a range in mag is provided
            q1, q2 = part.split('-')
            q1, q2 = int(q1), int(q2)
            quarters.append(q1)
            quarters.append(q2)
        else:
            # If only one mag-value is given select 1 mag around it
            q1 = int(part)
            quarters.append(q1)
    return quarters








def convertMagnitudeRange(dm):
    """
    Small function that takes a string of numbers (here of magnitudes)
    and split it up into readable float values used as real number
    ranges. If a single number is given, a selection of 1 mag around
    the imput int/float is returned as a magnitude range.
    Used in: PLATOnium/simulator-pic.py
    """
    magRange = []
    for part in dm.split(','):
        if '-' in part:
            # If a range in mag is provided
            m1, m2 = part.split('-')
            m1, m2 = float(m1), float(m2)
        else:
            # If only one mag-value is given select 1 mag around it
            m1 = float(part)-0.5
            m2 = float(part)+0.5
        magRange.append(m1)
        magRange.append(m2)
    return magRange







def getStarsWithinCameraGroup(camGroup, raPF, decPF, ra, dec, sizeSubfield=6):
    """
    This function determines if a star is within the FOV of a specific
    PLATO camera group. 

    Parameters
    ----------
    camGroup : int [1, 2, 3, 4]
        N-CAM camera group ID
    raPF : float
        Right acsension of pointing field [deg]
    decPF : float
        Declination of pointing field [deg]
    ra : list, array
        Right ascension of stars to be checked against [deg]
    dec : list, tuple, array
        Declination of stars to be checked against [deg]

    Return
    ------
    dex : list
       Indices booleans with True being all stars within camera group FOV
    distanceOA : list
       Distance of each from the camera's optical axis [mm]
    """

    # Setup for simulation object

    inputFile = os.getenv("PLATO_PROJECT_HOME") + "/inputfiles/inputfile.yaml"
    sim = Simulation(None, inputFile)

    # Telescope config

    quarter = 1

    raPlatformDeg  = sim["ObservingParameters/RApointing"]  = raPF   # [deg]
    decPlatformDeg = sim["ObservingParameters/DecPointing"] = decPF  # [deg]

    raPlatformRad  = np.deg2rad(raPlatformDeg)   # [rad]
    decPlatformRad = np.deg2rad(decPlatformDeg)  # [rad]

    focalLength      = float(sim["Camera/FocalLength/ConstantValue"]) * 1000.0  # [m] -> [mm]
    focalPlaneAngle  = np.deg2rad(float(sim["Camera/FocalPlaneOrientation/ConstantValue"]))

    solarPanelOrientation = sim["Platform/SolarPanelOrientation"] = math.fmod(quarter * 90., 360.) -6
    solarPanelOrientation = np.deg2rad(float(solarPanelOrientation))

    raTargetsRad  = np.deg2rad(ra)   # [rad]
    decTargetsRad = np.deg2rad(dec)  # [rad]

    # Loop over each star for this cam-group

    sim["Telescope/GroupID"] = camGroup
    azimuthTelescope = np.deg2rad(sim["CameraGroups/AzimuthAngle"][camGroup-1])
    tiltTelescope    = np.deg2rad(sim["CameraGroups/TiltAngle"][camGroup-1])

    dexGroup   = np.zeros(len(ra), dtype=bool)
    distanceOA = np.zeros(len(ra))

    for i in range(len(ra)):

        subfieldIsOnCCD = sim.setSubfieldAroundCoordinates(raTargetsRad[i], decTargetsRad[i],
                                                           sizeSubfield, sizeSubfield, normal=True)
        if subfieldIsOnCCD:

            xFP, yFP = skyToFocalPlaneCoordinates(raTargetsRad[i], decTargetsRad[i],
                                                  raPlatformRad, decPlatformRad,
                                                  solarPanelOrientation,
                                                  tiltTelescope, azimuthTelescope,
                                                  focalPlaneAngle, focalLength)

            distanceOA[i] = np.rad2deg(gnomonicRadialDistanceFromOpticalAxis(xFP, yFP,
                                                                             focalLength))

            if distanceOA[i] < 18.2: #sim['CCD/RelativeTransmissivity/RadiusFOV']:
                dexGroup[i] = True
            else:
                dexGroup[i] = False

        # Compile to bash
        compilation(i, len(ra), 'Group {}'.format(camGroup))
    print; print('')

    return dexGroup, distanceOA









def getInterPixelPositions(raPF, decPF, ra, dec):

    # Setup for simulation object

    inputFile = os.getenv("PLATO_PROJECT_HOME") + "/inputfiles/inputfile.yaml"
    sim = Simulation(None, inputFile)

    # Telescope config

    quarter = 1

    raPlatformDeg  = sim["ObservingParameters/RApointing"]  = raPF   # [deg]
    decPlatformDeg = sim["ObservingParameters/DecPointing"] = decPF  # [deg]

    raPlatformRad  = np.deg2rad(raPlatformDeg)   # [rad]
    decPlatformRad = np.deg2rad(decPlatformDeg)  # [rad]

    focalLength      = float(sim["Camera/FocalLength/ConstantValue"]) * 1000.0  # [m] -> [mm]
    focalPlaneAngle  = np.deg2rad(float(sim["Camera/FocalPlaneOrientation/ConstantValue"]))

    solarPanelOrientation = sim["Platform/SolarPanelOrientation"] = math.fmod(quarter * 90., 360.) -6
    solarPanelOrientationRad = np.deg2rad(float(solarPanelOrientation))

    raTargetsRad  = np.deg2rad(ra)   # [rad]
    decTargetsRad = np.deg2rad(dec)  # [rad]

    # Loop over each star for this cam-group

    camGroup = 1
    sim["Telescope/GroupID"] = camGroup
    azimuthTelescopeRad = np.deg2rad(sim["CameraGroups/AzimuthAngle"][camGroup-1])
    tiltTelescopeRad    = np.deg2rad(sim["CameraGroups/TiltAngle"][camGroup-1])

    # CCD properties

    pixelSize = float(sim["CCD/PixelSize"])
    
    ccdCode = np.zeros(len(ra))
    xCCD    = np.zeros(len(ra))
    yCCD    = np.zeros(len(ra))

    for i in range(len(ra)):

        subfieldIsOnCCD = sim.setSubfieldAroundCoordinates(raTargetsRad[i], decTargetsRad[i],
                                                           6, 6, normal=True)
        if subfieldIsOnCCD:
    
            # Fetch CCD code and pixel coordinates (account for field distortion in included)
            
            if sim["Camera/IncludeFieldDistortion"]:
                includeFieldDistortion = sim["Camera/IncludeFieldDistortion"]
                if sim["Camera/FieldDistortion/Type"] == 'FromFile':
                    mappedDistortion = True
                    distortionCoefficients = sim["Camera/FieldDistortion/CoefficientsFromFile"]
                else:
                    mappedDistortion = False
                    distortionCoefficients = sim["Camera/FieldDistortion/ConstantCoefficients"]
            else:
                includeFieldDistortion = False
                mappedDistortion = False
                distortionCoefficients = False
                
            ccdCode[i], xCCD[i], yCCD[i] = getCCDandPixelCoordinates(raTargetsRad[i], decTargetsRad[i],
                                                                     raPlatformRad, decPlatformRad,
                                                                     solarPanelOrientationRad,
                                                                     tiltTelescopeRad, azimuthTelescopeRad,
                                                                     focalPlaneAngle, focalLength, pixelSize,
                                                                     includeFieldDistortion, normal=True,
                                                                     mappedDistortion=mappedDistortion,
                                                                     distortionCoefficients=distortionCoefficients)

        # Compile to bash
        compilation(i, len(ra), 'Computing pixel positions')
    print; print('')
    
    return ccdCode, xCCD, yCCD








# def picOfDestiny(distribution, prange):
#     """
#     This function randomly picks a value from any gievn distribution and returns it.
#     The distribution must consist of values between 0 and 1 with its peak at 1. This function picks a
#     random value from the allowed range and then uses a distribution to get a P number
#     between 0 and 1, it then rolls a dice and chekcs wheter the dice roll is under the
#     P number. If it is, then the picked value is returned. This ensures a recration of
#     the distribution shape over thousands of picks
#     """

#     pick = random.random()*(prange[1]-prange[0]) + prange[0]
#     p = distribution(pick)
#     roll = random.random()
#     if roll < p:
#         return pick
#     else:
#         return distribution_pick(distribution, range)
