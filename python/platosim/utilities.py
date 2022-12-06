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

# PlatoSim
import platosim.referenceFrames as rf

#==============================================================#
#                           FUNCTIONS                          #
#==============================================================#


def errorcode(API, message):

    """Function to colour code error messages within a code.

    Parameters
    ----------
    API : str
       Which API to use: [software, module, message, warning, error]
    message : str
       Message to add to API

    Return
    ------
    Error message written to bash
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





def tqdm_bar_format():

    """Code snippet to set default
    """

    bar_format = "{l_bar}{bar:50}{r_bar}{bar:-50b}"
    return bar_format

        


        
def compilation(i, i_max, text=''):

    """Custum function to print out a compilation-time-bar in the terminal.

    Parameters
    ----------
    i : int
        Index of the current running job-loop
    i_max : int
        Index of the last running job-loop
    text : str
        Optional text written next to compilation-bar

    Return
    ------
        No parameters only a nice compilation bar to bash

    Example
    -------
    for i in range(10):
        compilation(i, len(ra), 'Computing pixel positions')
    print; print('')
    """

    # Running percentage calculation

    percent = (i + 1) / (i_max * 1.0) * 100

    # We here divide by 2 as the length of the bar is only 50 characters:

    bar = "[" + "-" * int(percent/2) + '>' + " " * (50-int(percent/2))+"] {}% {}"\
        .format(int(percent), text)

    # Print and clean-up

    sys.stdout.write(u"\u001b[1000D" + bar)
    sys.stdout.flush()






def medianAbsoluteDeviation(array):
    """
    Calculate the Median Abosolute Deviation (MAD) of an array.
    """
    return np.sum( np.abs(array - np.median(array)) ) / len(np.ravel(array))



def rootMeanSquare(array):
    """
    Calculate the Root Mean Square (RMS) of an array.
    """
    return np.sqrt( np.mean(array**2) )


    


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

    Parameters
    ----------
    filtType : string
        Filter is either: median or mean
    signalIn : ndarray
       Signal needed for processing
    numBox : int
        Integer used as car-box size of 1 hour: 3600s/25s = 144.

    Return
    ------
    signalOut : ndarray
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
    
    """Coversion from Johnson-Cousin V magnitude to the PLATO passband.
    
    This filtersion is from Marchiori et al. (2019), Eq. 5 and 6, and is
    extracted using synthetic stellar spectra from the POLLOX database.
    NOTE valid for Teff = 4000-15000 K (hence not for M-dwarfs).

    Parameters
    ----------
    V : float, narray
        Johnson-Cousin magnitude of star(s).
    Teff : float narray
        Effective temperature of star(s).

    Return
    ------
    P : float, narray
        The PLATO passband magnitude of star(s).
    """

    # The actual filtersion equation

    c = [1.184e-12, 4.526e-8, 5.805e-4, 2.449]     # Machiori et al. (2019)
    #c = [2.366e-12, 8.126e-08, -0.0009279, 3.499] # Fabio Fialho et al. in prep
    P  = c[0]*Teff**3 - c[1]*Teff**2 + c[2]*Teff - c[3] + V

    return P







def getPhotonNoiseLimitNSR(P, Ncam=1, Ntra=1, tdur=3600, camType='N'):

    """NSR estimate in the photon noise limit of bright stars. 

    The stellar flux are calculated from the PLATO passband found by 
    Marchiori et al. (2019).
    
    NOTE: only valid for very bright stars (P < 11).

    Parameters
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

    Return
    ------
    NSR : float, narray
        NSR only valid for the photon noise limit.
    """

    # Choose cycle and exposure time [s] for either the normal (N) or fast (F) cameras

    if camType == 'N':
        texp = 21.
        tcyc = 25.
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

    """Function to add a column to an exisiting pandas data frame.
    """
    
    df[name] = newCol
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    return df[cols]






def convertQuarterRange(dQ):

    """Function to sort a quarter ranges.
    
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







def stellarFlux(Vmag, exposureTime, fluxm0=1.00238e8,
                throughputBandwidth=400, transmissionEfficiency=0.76, 
                lightCollectingArea=0.01131, quantumEfficiency=0.87):

    """
    PURPOSE: compute the stellar flux (electrons / exposure) given the instrumental characteristics

    INPUT: Vmag:                   Johnson V magnitude
           exposureTime:           Exposure time (without the readout) [s]
           fluxm0:                 Photon flux of a V=0 star (default SpT=G2V) [phot/s/m^2/nm]
           throughputBandwidth:    FWHM [nm]
           transmissionEfficiency: In [0,1]
           lightCollectingArea:    Of the telescope [m^2]
           quantumEfficiency:      In [0,1]

    OUTPUT: flux: [e-/exposure]
    """

    photonFlux = (fluxm0 * throughputBandwidth * transmissionEfficiency *
                  lightCollectingArea * pow(10.0, -0.4 * Vmag) * exposureTime)
    electronFlux = photonFlux * quantumEfficiency

    return electronFlux







def convertMagnitudeRange(dm):

    """Function to sort magnitudes ranges.

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








def imageNorm(inputArray, norm="linear", sigma=2, scale_min=None, scale_max=None):
    """
    Performs custom scaling of the input numpy array.

    @type inputArray: np array
    @param inputArray: image data array
    @type scale_min: float
    @param scale_min: minimum data value
    @type scale_max: float
    @param scale_max: maximum data value
    @rtype: np array
    @return: image data array
    """
    # Input image array
    
    image = np.array(inputArray, copy=True)

    # Default scaling is 2 sigma

    if scale_min is None:
        scale_min = image.mean() - sigma*image.std()
    if scale_max is None:
        scale_max = image.mean() + sigma*image.std()

    # Clip data
    
    image = image.clip(min=scale_min, max=scale_max)
    
    # Select normalization method
    
    if norm == "linear":
        image   = (image - scale_min) / (scale_max - scale_min)
        indices = np.where(image < 0)
        image[indices] = 0.0
        indices = np.where(image > 1)
        image[indices] = 1.0

    elif norm == "log":
        factor = np.log10(scale_max - scale_min)
        indices0 = np.where(image < scale_min)
        indices1 = np.where((image >= scale_min) & (image <= scale_max))
        indices2 = np.where(image > scale_max)
        image[indices0] = 0.0
        image[indices2] = 1.0
        image[indices1] = np.log10(image[indices1]) / factor
        
    elif norm == "sqrt":
        image = image - scale_min
        indices = np.where(image < 0)
        image[indices] = 0.0
        image = np.sqrt(image)
        image = image / math.sqrt(scale_max - scale_min)
        image = np.sqrt(image)
        image = image / np.sqrt(scale_max - scale_min)

    elif norm == "asinh":
        non_linear = 2.0
        factor = np.arcsinh((scale_max - scale_min) / non_linear)
        indices0 = np.where(image < scale_min)
        indices1 = np.where((image >= scale_min) & (image <= scale_max))
        indices2 = np.where(image > scale_max)
        image[indices0] = 0.0
        image[indices2] = 1.0
        image[indices1] = np.arcsinh( (image[indices1] - scale_min) / non_linear) / factor

    #else:
    #    errorcode("error", "Not valid normalization method!")

    # Finito!
        
    return image









def moveColorbarExponent(x_offs=0, y_offs=1, dig=0, side='left', omit_last=False):

    """Move scientific notation exponent from top to the side.
    
    Additionally, one can set the number of digits after the comma
    for the y-ticks, hence if it should state 1, 1.0, 1.00 and so forth.

    Parameters
    ----------
    offs : float, optional; <0>
        Horizontal movement additional to default.
    dig : int, optional; <0>
        Number of decimals after the comma.
    side : string, optional; {<'left'>, 'right'}
        To choose the side of the y-axis notation.
    omit_last : bool, optional; <False>
        If True, the top y-axis-label is omitted.

    Returns
    -------
    locs : list
        List of y-tick locations.

    Note
    ----
    This is kind of a non-satisfying hack, which should be handled more
    properly. But it works. Functions to look at for a better implementation:
    ax.ticklabel_format
    ax.yaxis.major.formatter.set_offset_string
    """

    # Get the ticks
    locs, _ = plt.yticks()

    # Put the last entry into a string, ensuring it is in scientific notation
    # E.g: 123456789 => '1.235e+08'
    llocs = '%.3e' % locs[-1]

    # Get the magnitude, hence the number after the 'e'
    # E.g: '1.235e+08' => 8
    yoff = int(str(llocs).split('e')[1])

    # If omit_last, remove last entry
    if omit_last:
        slocs = locs[:-1]
    else:
        slocs = locs

    # Set ticks to the requested precision
    form = r'$%.'+str(dig)+'f$'
    plt.yticks(locs, list(map(lambda x: form % x, slocs/(10**yoff))))

    # Define offset depending on the side
    if side == 'left':
        x_offs = -.18 - x_offs # Default left: -0.18
    elif side == 'right':
        x_offs = 1 + x_offs    # Default right: 1.0
        
    # Plot the exponent
    plt.text(x_offs, y_offs, r'$\times10^{%i}$' % yoff, transform =
            plt.gca().transAxes, verticalalignment='top')

    # Return the locs
    return locs








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



