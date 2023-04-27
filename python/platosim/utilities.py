#!/usr/bin/env python3

"""
This python module contains all general utilities that are commonly used
by the different codes within PlatoSim and PLATOnium.

NOTE: these utilities needs the Poetry install!
"""

# Standard
import os
import sys
import glob
import math
import ftplib
import inspect
import fnmatch

# Extra
import h5py
import pathlib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pylab import MaxNLocator
from prettytable import PrettyTable
from scipy.ndimage import median_filter

# PlatoSim
import platosim.referenceFrames as rf


#--------------------------------------------------------------#
#                        BASH FUNCTIONS                        #
#--------------------------------------------------------------#


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
    from colorama import Fore, Style
    
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
        sys.exit()





def fileMatch(fileList, stringList):
    """
    fileMatch(fileList, stringList)
    
    returns the list of files from fileList in which all strings from stringList 
    can be found, in any order Case sensitive.
    """

    while (len(stringList)>1):
        fileList = fileMatch(fileList,[stringList.pop()])
    return [file for file in fileList if fnmatch.fnmatch(file,'*'+stringList[0]+"*")]





def fileSelect(stringList, location="./", listOrder=0):

    """fileSelect(stringList, location="./", listOrder=0)

    Returns a list of all files in 'location' with name matching every string in the list
    If listOrder is True, the ordering is forced to be identical to the one in stringList
    """
    
    allfiles = os.listdir(location)
    if listOrder:
        pattern = "*"
        i = 0
        while i < len(stringList): 
            pattern += stringList[i] + "*"
            i += 1
        return [file for file in allfiles if fnmatch.fnmatch(file,pattern)]
    else:
        return fileMatch(allfiles, stringList)




    
def getFunctions(script):

    """Fetch names of class functions.
    """
    names = inspect.getmembers(script, inspect.isfunction)
    funcs = [item[0] for item in names]
    t = PrettyTable()
    t.add_column(f"{script.__name__} functions", funcs)
    return t
    




def tqdmBar():

    """Add-on function to be used in a for statement with tqdm library.

    Example:
    >>> from tqdm import tqdm
    >>> for i in tqdm(range(<number of loops>), bar_format=ut.tqdmBar()):
            <loop over something>
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

    bar = "[" + "-" * int(percent / 2) + '>' + " " * (50 - int(percent / 2)) + "] {}% {}" \
        .format(int(percent), text)

    # Print and clean-up

    sys.stdout.write(u"\u001b[1000D" + bar)
    sys.stdout.flush()


    


def downloadFromFTP(filename, outputDir, server='plato'):

    """Function to download file from KUL FTP.
    https://stackoverflow.com/questions/67300881/how-do-i-keep-a-ftp-connection-alive
    """

    # Assume that no suffix means a folder of data
    # If true then download folder and it entire content
    # If flase simply download the requested file
    
    ftp_filename = pathlib.Path(filename)

    if ftp_filename.suffix in ('.zip', '.npy', '.ftr', '.hdf5', '.h5'):
        ftp_subpath = pathlib.Path(filename).parents[0]
        permission  = False
    else:
        ftp_subpath = pathlib.Path(filename)
        permission  = True
        
    # Also if file on FTP is within a folder, create folder locally
        
    outputDir = outputDir / ftp_subpath
    outputDir.mkdir(parents=True, exist_ok=True)
        
    # Login to server
    # For plato: Download a single file
    # For platodata: Download all files in a folder

    ftp = ftplib.FTP('ftp.ster.kuleuven.be')
    
    if server == 'plato':
        # Single file download
        ftp.login(user=server, passwd='miSotalP')
        ftp.cwd(f'{ftp_subpath}')
        files = [filename]
        #ftp = 'ftp://plato:miSotalP@ftp.ster.kuleuven.be'
    elif server == 'platodata':
        ftp.login(user=server, passwd='i9Pidw1bXIFShGYb0jI8')
        ftp.cwd(f'PLATOSIM/{ftp_subpath}')
        files = ftp.nlst()[2:]
        #ftp = 'ftp://platodata:i9Pidw1bXIFShGYb0jI8@ftp.ster.kuleuven.be/PLATOSIM'
    else:
        errorcode('error', f'Server name {server} is not valid!')
            
    # Fetch all the files
        
    for filename in files:

        # Only try to save file if is doesn't exists
        
        local_file = pathlib.Path(outputDir) / filename

        if not local_file.is_file():
            ftp_file   = open(local_file, 'wb')
            ftp.retrbinary(f'RETR {filename}', ftp_file.write)
            ftp_file.close()
            
            # Give read and write rights to this
            if permission: local_file.chmod(777)

    # Close connection
    
    ftp.quit()

    

            
#--------------------------------------------------------------#
#                      PANDAS OPERATIONS                       #
#--------------------------------------------------------------#


def pdAddColumn(df, newCol, name):

    """Add a column to an exisiting pandas data frame as first entry.
    """

    df[name] = newCol
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    return df[cols]





#--------------------------------------------------------------#
#                       NUMPY OPERATIONS                       #
#--------------------------------------------------------------#

    
def findNearestIndex(array, value):

    """Find the nearest value within an numpy array.
    """

    return (np.abs(np.asarray(array) - value)).argmin()





def medianAbsoluteDeviation(array):

    """Calculate the Median Abosolute Deviation (MAD) of an array.
    """
    
    return np.sum(np.abs(array - np.median(array))) / len(np.ravel(array))





def rootMeanSquare(array):

    """Calculate the Root Mean Square (RMS) of an array.
    """
    
    return np.sqrt(np.mean(array ** 2))





def normalize(signal, factor=1e6, length=-1):

    """Normalize a signal with the option to factorize it.

    Parameters
    ----------
    signal : narray
        Signal array that needs to be turned into relative signal.
    factor : int, float
        Factor to scale signal. Default is in ppm.
    length : int, float
        Option to normalize signal using only a smaller initial signal segment.
        Default is to use the entire signal sequence.

    Return
    ------
    relative_signal : narray
        Normalized relative signal is returned. Default unit in [ppm].
    """

    relative_signal = (signal / np.nanmean(signal[:int(length)]) - 1) * factor

    return relative_signal





def imageNorm(inputArray, norm="linear", sigma=2, scale_min=None, scale_max=None):

    """Performs custom scaling of the input numpy array.

    Parameters
    ----------
    inputArray : ndarray
        Input image array to normalize.
    norm : str
        Normalization method. 
        Options: ['linear', 'log', 'sqrt', 'asinh'] 
    sigma : float
        Scaling factor corresponding to the std of the image.
    scale_min : float
        Minimum data value.
    scale_max : float
        Maximum data value.
    
    Return
    ------
    image : ndarray
        Normalized image array.
    """
    
    # Input image array

    image = np.array(inputArray, copy=True)

    # Default scaling is 2 sigma

    if scale_min is None:
        scale_min = image.mean() - sigma * image.std()
    if scale_max is None:
        scale_max = image.mean() + sigma * image.std()

    # Clip data

    image = image.clip(min=scale_min, max=scale_max)

    # Select normalization method

    if norm == "linear":
        image = (image - scale_min) / (scale_max - scale_min)
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
        image[indices1] = np.arcsinh((image[indices1] - scale_min) / non_linear) / factor

    # That's it!

    return image




#--------------------------------------------------------------#
#                        PLATO SPECIFIC                        #
#--------------------------------------------------------------#


def stellarFlux(Vmag, exposureTime, fluxm0=1.00238e8,
                throughputBandwidth=400, transmissionEfficiency=0.76,
                lightCollectingArea=0.01131, quantumEfficiency=0.87):

    """Compute the stellar flux given the instrumental characteristics.

    Parameters
    ----------
    Vmag : float
        Johnson-Cousin V magnitude
    exposureTime : float
        Exposure time (without the readout) [s]
    fluxm0 : float
        Photon flux of a V=0 G2V star [phot/s/m^2/nm]
    throughputBandwidth : float 
        Throughput FWHM value in the passband [nm]
    transmissionEfficiency : 
        Transmission efficiency in the passband [0, 1]
    lightCollectingArea : float
        The aperture of the cameras [m^2]
    quantumEfficiency : float
        Quantum efficiency of the detector [0, 1]

    Return
    ------
    flux : float
        Instrumental stellar flux [e-/exposure]
    """

    photonFlux = (fluxm0 * throughputBandwidth * transmissionEfficiency *
                  lightCollectingArea * pow(10.0, -0.4 * Vmag) * exposureTime)
    electronFlux = photonFlux * quantumEfficiency

    return electronFlux





def fromMagToRelativeFlux(mag, norm=1e6):

    """Convert magnitude to relative flux

    Parameters
    ----------
    mag : float
        Input magnitude
    norm : float
        Normalisation contant for relative flux

    Return
    ------
    flux : ndarray
        Relative flux scaled after the normalisation constant.
    """

    return (10**(-0.4*mag) - 1) * norm






def passbandConversionV2P(mag, Teff, inverse=False, method='fialho'):

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

    # Bolometric scaling relation

    if method == 'fialho':
        # Fialho et al. (in prep.)
        c   = [-2.366e-12, 8.126e-8, -9.279e-4, 3.499]
    elif method == 'marchiori':
        # Machiori et al. (2019)
        c = [-1.184e-12, 4.526e-8, -5.805e-4, 2.449]

    bol = c[0]*Teff**3 + c[1]*Teff**2 + c[2]*Teff + c[3]

    # From V to P (or P to V if inverse it True)

    if inverse:
        return mag + bol
    else:
        return mag - bol





def getJitterNoiseLimitNSR(rms, tdur=3600, level='instrument', camType='normal'):

    """NSR estimate of the jitter noise component.

    Parameters
    ----------
    jitterASD : float
        Amplitude Spectral Density of the jitter [ppm/sqrt(muHz)] :
        If level = 'camera'     : At the cycle frequency of the cameras
        If level = 'instrument' : Over the duration of all exposures
    tdur : float, narray
        Time duration over which the NSR is estimated. E.g., 3600s for 1h precision.
    camType : str
        Either the normal (N) or fast (F) cameras. Default is normal.

    Return
    ------
    NSR : float
        NSR only valid for the photon noise limit.
    """

    # Amplitude Spectral Density of the jitter [ppm/sqrt(muHz)]
    
    jitterASD = rms / np.sqrt(1e-6)

    # Choose cycle and exposure time [s] for either the normal (N) or fast (F) cameras

    if camType == 'normal':
        tcyc = 25.
    elif camType == 'fast':
        tcyc = 2.5

    # Number of images to average over

    nimg = int(tdur/tcyc)

    # Calculate the jitter noise

    if level == 'camera':
        jitterNoise = jitterASD * np.sqrt(1 / tcyc)
    elif level == 'instrument':
        jitterNoise = jitterASD * np.sqrt(1 / (tcyc*nimg) ) 
    else:
        errorcode('error', 'No such "level" entry!')
        
    # Return

    return jitterNoise




        
def getPhotonNoiseLimitNSR(mag, passband='P', camType='normal', ncam=1, ntra=1, tdur=3600):

    """NSR estimate in the photon noise limit of bright stars.

    The stellar flux are calculated from the PLATO passband found by 
    Marchiori et al. (2019).

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

    # Choose cycle and exposure time [s] for either the normal or fast cameras

    if camType == 'normal':
        texp = 21.
        tcyc = 25.
        gain = 0.0222 * 2.14   # [ADU/e-]
    else:
        texp = 2.1
        tcyc = 2.5
        gain = 0.05   # TODO: update gain values for F-CAMs

    # Flux of stars [e-/s]
    
    if passband == "V":
        f0 = 1.00179e8
        f = 10**(-0.4 * mag) * f0
        
    elif passband == 'P':
        # The P passband zero-point
        if camType == 'normal':
            zp   = 20.77
        if camType == 'fastblue':
            zp = 20.18
        if camType == 'fastred':
            zp = 19.81
        # Calculate flux
        f0 = 0.7324478224428527e8
        f = 10**(-0.4 * (mag - zp)) * f0
    else:
        errorcode('error', f'Wrong {camType} name!')

    # Observed total flux [ADU/exp]

    F = f * tcyc * gain

    # SNR from pure photon noise and NSR from uncorrelated noise.
    # Gaussian statistic gives sigma --> sigma/sqrt(N)
    
    NSR = 1e6 / np.sqrt(F * ncam * ntra * tdur)

    return NSR





def getBackgroundNoiseLimitNSR(mag, passband='P', camType='normal', tdur=3600):

    """NSR estimate in the photon noise limit of bright stars.

    The stellar flux are calculated from the PLATO passband found by 
    Marchiori et al. (2019).

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

    # Choose cycle and exposure time [s] for either the normal or fast cameras

    if camType == 'normal':
        gain = 1/(0.0222 * 2.14)   # [ADU/e-/pixel]
    else:
        gain = 0.05   # TODO: update gain values for F-CAMs

    if passband == "V":
        f0 = 1.00179e8
    elif passband == 'P':
        f0 = 0.7324478224428527e8
    else:
        errorcode('error', f'Wrong {camType} name!')
        
    # Calculate noise and signal
    gain = 25
    bg   = 60  # [e-/s/pixel]
    mask = 20  # [pixel]
    throughput   = 0.8134999994206865
    transmission = 0.4822896122932434
    
    noise  = gain * bg * tdur * mask * throughput #* transmission
    signal = np.sqrt(10**(-0.4 * mag) * f0 * tdur)**1.8
    return noise / signal



#--------------------------------------------------------------#
#                       PLATOnium FUNCTIONS                    #
#--------------------------------------------------------------#


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

    # That's it!
            
    return quarters





def convertMagnitudeRange(dm):

    """Function to sort magnitudes ranges.

    Small function that takes a string of numbers (here of magnitudes)
    and split it up into readable float values used as real number
    ranges. If a single number is given, a selection of 1 mag around
    the imput int/float is returned as a magnitude range.
    
    Used in: PLATOnium/simulator-pic.py

    Parameters
    ----------
    dm : str
        Magnitude value or range (e.g. '10.5' or '10.0-11.5')

    Return
    ------
    magRange : range()
        Corresponding magnitude range using the range() object.
    """
    
    magRange = []
    for part in dm.split(','):
        
        if '-' in part:

            # If a range in mag is provided

            m1, m2 = part.split('-')
            m1, m2 = float(m1), float(m2)
            
        else:
            
            # If only one mag-value is given select 1 mag around it
            
            m1 = float(part) - 0.5
            m2 = float(part) + 0.5
            
        magRange.append(m1)
        magRange.append(m2)

    # That's it!
        
    return magRange        




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

