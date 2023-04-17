# -*- coding: utf-8 -*-

"""
TODO this file is really not needed since its functionalities has been 
implemented into simfile.py. Maybe check which method is best.
"""

import numpy as np
import matplotlib.cm as cm
from matplotlib import pyplot as plt 
    
import platosim.utilities import ut



class Formatter:
    """
    This formatter class is used to make matplotlib show the pixel coordinates and intensity 
    of an image shown by imshow()
    """

    def __init__(self, im):
        self.im = im
    def __call__(self, x, y):
        z = self.im[int(y), int(x)]
        return 'x={:.01f}, y={:.01f}, z={:.01f}'.format(x, y, z)












def showSim(simfile, n=0, type="Image", dpi=0, figname=None,verbose=True, **kwargs):
    
    """
    SYNOPSIS
        showSimData(simfile, n=0, type="Image", dpi=0,figname=None, verbose=1, **kwargs)

    INPUTS
        simfile  : platosim file opened with h5py.File
        n = 0    : image number to display
        type     : first 2 characters are used to branch to 
                   im - pixelImage, 
                   su - subPixelImage
                   sm - smearingMap
                   bi - biasMap
                   fp - prnu -- FF pixel
                   fs - irnu -- FF subpixel

        dpi      : for subpixelImage and IRNU, controls the image resolution (see output)
        figname  : change the default matplotlib.figure name
        verbose  : False : silent. True : prints clim (except for Flatfield)
        **kwargs : any additional matplotlib.pyplot.imshow keyword arguments
   
    OUTPUT
        matplotlib.imshow with following defaults
                    cmap          : hot
                    interpolation : nearest
                    origin        : lower
                    dpi for subpixel data : 96

    EXAMPLES:
        >>> from showSim import showSim
        >>> import h5py
        >>> myFile = h5py.File(path, "r")
        >>> showSim(myFile, 3, type="im")
        >>> showSim(myFile, type="pr", figname="Flatfield")
    """
    from platosim.h5 import h5get

    # Extract the data from the HDF5 file

    image = getSim(simfile,type=type,n=n)
    
    # If the user didn't specify the figure name, create our own figure name that is based
    # on the requested data product. 

    firstTwoLettersOfDataProduct = str(type).lower()[:2]
    if not figname:
       dataProduct = {"im":'image', "su":'subPixelImage', "bi":"biasMap",\
       "sm":"smearingMap", "fp":"PRNU", "fs":"IRNU", "pp":"rebinnedPSFpixel",\
       "ps":"rebinnedPSFsubPixel", "pr":"rotatedPSF"}
       figname = dataProduct[firstTwoLettersOfDataProduct]+str(n).zfill(6)
    
    # Compute the median and the deviations of the image, which will be used to set the intensity range
    # in the imshow plot.
    
    med    = np.median(image)
    meddev = np.max([1., ut.mad(image)])
    stddev = np.max([1., np.std(image)])
    
    # Set the default imshow arguments

    kwargs.setdefault('interpolation', 'nearest')
    kwargs.setdefault("origin", "lower")
    if firstTwoLettersOfDataProduct not in ["pp","pr","ps"]: 
        kwargs.setdefault('cmap', cm.hot)
    if firstTwoLettersOfDataProduct not in ["fp", "fs","pp","pr","ps"]: 
        kwargs.setdefault("clim",[int(med-meddev),int(med+3.*stddev)])
        if verbose: print("clim", kwargs["clim"])
    
    # Show the plot
    """
    # JORIS' VERSION : ALWAYS CREATES A NEW PLOT
    if firstTwoLettersOfDataProduct in ["su","ir"]:
        if not dpi: dpi = 96
        fig, axis = plt.subplots(dpi=dpi)
    else:
        fig, axis = plt.subplots()
    """
    # PIERRE'S VERSION : UNLESS SPECIFIED VIA FIGNAME, REUSES EXISTING PLOT     
    if firstTwoLettersOfDataProduct in ["su","ir"]:
        if not dpi: dpi = 96
        fig = plt.figure(figname,dpi=dpi)
        axis = fig.add_subplot(111)
    else:
        fig = plt.figure(figname)
        axis = fig.add_subplot(111)

    
    zeroPointRow = h5get(simfile,"ZeroPointRow",verbose=False)
    zeroPointColumn = h5get(simfile,"ZeroPointColumn",verbose=False)
    
    axis.set_title(figname + " - LL Corner: [{0},{1}]".format(zeroPointRow,zeroPointColumn))
    axis.imshow(image, **kwargs)
    plt.show()

    # That's it!

    return axis








def getSim(simfile, n=0, type="Image"):
    
    """
    SYNOPSIS
        getSim(simfile, n=0, type="Image")
        From the HDF5 file, extract from the image data of the specified type in a numpy array.

    INPUTS
        simfile  : platosim file opened with h5py.File
        n = 0    : image number to display
        type     : first 2 characters are used to branch to 
                   im - pixelImage, 
                   su - subPixelImage
                   sm - smearingMap
                   br - biasMapRight
                   bl - biasMapLeft
                   fp - prnu -- FF pixel
                   fs - irnu -- FF subpixel
                   pp - rebinnedPSFpixel
                   ps - rebinnedPSFsubPixel
                   pr - rotatedPSF
        
    OUTPUT
        Numpy arrays with the requested data product according to input parameters 'n' and 'type'
        Flatfield PRNU and IRNU are unique (=> n is ignored in these cases)
    """

    # The data arrays are stored within certain groups in the hdf5 file. The following is to construct
    # the group names and the data array name within that group

    firstTwoLettersOfDataProduct = str(type).lower()[:2]
    hdf5GroupNameMapping = {"im":'Images', "su":'SubPixelImages', "br":"BiasMapsRight", "bl":"BiasMapsLeft",\
    "sm":"SmearingMaps", "fp":"Flatfield", "fs":"Flatfield","pp":"PSF","ps":"PSF","pr":"PSF"}
    hdf5DataNameMapping  = {"im":'image', "su":'subPixelImage', "br":"biasMap", "bl":"biasMap",\
    "sm":"smearingMap", "fp":"PRNU", "fs":"IRNU", "pp":"rebinnedPSFpixel", "ps":"rebinnedPSFsubPixel", "pr":"rotatedPSF"}
    
    groupName = hdf5GroupNameMapping[firstTwoLettersOfDataProduct]

    # Extract the data into numpy arrays

    if firstTwoLettersOfDataProduct in ["fp","fs","pp","ps","pr"]:
        dataName = hdf5DataNameMapping[firstTwoLettersOfDataProduct]
        return simfile[groupName][dataName]
    elif firstTwoLettersOfDataProduct in ["im","su","bl", "br","sm"]:
        dataName = hdf5DataNameMapping[firstTwoLettersOfDataProduct] + str(n).zfill(6)
        return simfile[groupName][dataName]
    else:
        raise Exception("Unrecognised type. See accepted types in getSimData.__doc__. Type: {}".format(type))

