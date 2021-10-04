#!/usr/bin/env python3

from numpy import *
import numpy as np
from scipy.ndimage import median_filter

from matplotlib import pyplot as plt
from matplotlib import patches
from matplotlib.path import Path
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

import astropy.units as u
from astropy.coordinates import SkyCoord

from platosim.photometryfile import PhotometricFile
from platosim.referenceFrames import *
from platosim.utilities import *

# Top level Matplotlib settings to ease the writing

fs = 15    # Font size
lw = 0.3   # Line width
ms = 0.5   # Scatter plot size

# Constants

rad2arcsec = 648000 / np.pi
day2sec    = 86400.

#==============================================================#
#                         GRAPHICAL TOOLS                      #
#==============================================================#

def axes_minmax(x=None, y=None, pt=0.02):
    """
    This is a small utility to automatically scale the axes of a plot
    using a default spacing in percentage (pt). The user need only to
    specify which axis a min and max limit should be returned from.
    """
    if x is not None:
        axmin = x[0]  - (x[-1]-x[0])*pt
        axmax = x[-1] + (x[-1]-x[0])*pt
    if y is not None:
        axmin = np.min(y) - (np.max(y)-np.min(y))*pt
        axmax = np.max(y) + (np.max(y)-np.min(y))*pt
    return axmin, axmax


def axes_maskupdates(ax, time, maskupdates):
    """
    This is a small utility that takes an axes object, time points
    from a time series, and the mask-updates given in the same unit
    of time as the time points, and then plots vertical lines for
    every mask-update and quarter marks.
    """

    # Plot occurance of mask update

    #updates = np.arange(0, time[-1], maskupdate)
    for update in maskupdates:
        if update == 0:
            ax.axvline(x=update, c='k', linestyle=':', linewidth=1, label='Mask updates')
        else:
            ax.axvline(x=update, c='k', linestyle=':', linewidth=1)

    # Plot quarters

    quarters = np.arange(0, time[-1], 90)
    for Q in quarters:
        if Q == 0:
            ax.axvline(x=Q, c='darkgray', linestyle='-.', linewidth=1, label='Quarter marks')
        else:
            ax.axvline(x=Q, c='darkgray', linestyle='-.', linewidth=1)


#==============================================================#
#                      GRAPHICAL FUNCTIONS                     #
#==============================================================#


def drawCCDsInSkyMollweide(fig, raPlatform, decPlatform, solarPanelOrientation, tiltAngle, azimuthAngle, focalPlaneAngle, focalLength, pixelSize, normal=True):

    """
    PURPOSE: Project and plot the 4 CCDs of 1 camera on the sky

    INPUT: raPlatform:            right ascension of the platform pointing axis             [rad]
           decPlatform:           declination of the platform pointing axis                 [rad]
           solarPanelOrientation: (0,pi/2,pi,3pi/2) for quarters (Q1,Q2,Q3,Q4)              [rad]
           tiltAngle:             tilt angle of the telescope w.r.t. platform z-axis        [rad]
           azimuthAngle:          azimuth angle of the telescope on the platform            [rad]
           focalPlaneAngle:       angle between the Y_TL axis and the Y_FP axis: gamma_FP   [rad]
           focalLength:           focal length of the camera                                [mm]
           pixelSize:             pixel size                                                [micron]
           normal:                True for the normal camera configuration, False for the fast cameras

    OUTPUT: None

    TODO: - Does not work yet for the fast cams
          - Does not take distortion into account yet
    """

    # Select the proper CCD codes depending on whether we're dealing with the nominal or the fast cams

    if normal == True:
        ccdCodes = ['1', '2', '3', '4']
    else:
        ccdCodes = ['1F', '2F', '3F', '4F']


    # Set up the colors to be used to draw each CCD.
    # Different CCDs have different colors.

    color = {'1': 'b', '1F': 'b', '2': 'r', '2F': 'r', '3': 'g', '3F': 'g', '4': 'k', '4F': 'k'}

    # Set up the figure

    axes = fig.add_subplot(111, projection="mollweide")
    axes.grid(True)

    # Plot each of the 4 CCDs

    for ccdCode in ccdCodes:

        # Get the focal plane FP' coordinates of the CCD corners  [mm]

        cornersXmm, cornersYmm = computeCCDcornersInFocalPlane(ccdCode, pixelSize)

        # Compute the equatorial sky coordinates [rad] from the the focal plane FP' coordinates [mm] of the corners

        ra, dec = focalPlaneToSkyCoordinates(cornersXmm, cornersYmm, raPlatform, decPlatform, solarPanelOrientation,  \
                                             tiltAngle, azimuthAngle, focalPlaneAngle, focalLength)

        # Repeat the coordinates of the 1st corner, to plot a nice closed loop
        # Convert from radians to degrees

        ra  = append(ra, ra[0]) 
        dec = append(dec, dec[0]) 

        # The sky projection expects a longitude in [-pi, +pi] rather than [0, 2* pi]
        # Moreover, the longitude should be reversed so that East is to the left
        
        ra[ra>pi] -= 2*pi
        ra = -ra 

        axes.plot(ra, dec, c=color[ccdCode])

        # Overplot the row closest to the readout register with a thicker line

        axes.plot([ra[0], ra[1]], [dec[0], dec[1]], c=color[ccdCode], linewidth=3)


    # Change the tick labels so that they are 0->360, rather than -180->+180

    tickLabels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tickLabels = np.remainder(tickLabels+360, 360)
    axes.set_xticklabels(tickLabels)     

    # Add axis labels

    plt.xlabel("RA [deg]")
    plt.ylabel("Dec [deg]")
    plt.draw()

    # That's it

    return axes







def drawStarsInSkyMollweide(fig, ra, dec):

    """
    PURPOSE: Project and plot the stars with the given right ascension and declination on the sky

    INPUT: ra:      right ascension of the stars             [degrees]
           dec:     declination of the areA                 [degrees]

    OUTPUT: None
    """
    
    # Set up the figure

    axes = fig.add_subplot(111, projection="mollweide")
    axes.grid(True)

    raRadians = []
    decRadians = []
    
    for index in range(len(ra)):
        raRadians.append(-ra[index] * pi / 180.0)
        decRadians.append(dec[index] * pi / 180.0)

    axes.plot(raRadians, decRadians, 'ko')


    # Change the tick labels so that they are 0->360, rather than -180->+180

    tickLabels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tickLabels = np.remainder(tickLabels+360, 360)
    axes.set_xticklabels(tickLabels)     

    # Add axis labels

    plt.xlabel("RA [deg]")
    plt.ylabel("Dec [deg]")
    plt.draw()

    # That's it

    return axes











def drawCCDsInFocalPlane(pixelSize, plotCCDlabels=True, normal=True):

    """
    PURPOSE: Plot the 4 CCDs in the focal plane in the FP' reference frame.
             May serve as a background to overplot the projected stars on the focal plane

    INPUT: pixelSize: size of 1 pixel [micron]
           normal: True for the normal camera configuration, False for the fast cameras

    OUTPUT: None

    """

    # Select the proper CCD codes depending on whether we're dealing with the nominal or the fast cams
    
    if normal == True:
        ccdCodes = ['1', '2', '3', '4']
    else:
        ccdCodes = ['1F', '2F', '3F', '4F']


    # Set up the colors to be used to draw each CCD. 
    # Different CCDs have different colors.

    color = {'1': 'b', '1F': 'b', '2': 'r', '2F': 'r', '3': 'g', '3F': 'g', '4': 'k', '4F': 'k'}


    # Plot each of the 4 CCDs

    fig = plt.figure(figsize = (10,10))
    ax = fig.add_subplot(111)

    for ccdCode in ccdCodes:

        # Get the corner coordinates in the FP' plane

        cornersXmm, cornersYmm = computeCCDcornersInFocalPlane(ccdCode, pixelSize)

        # Repeat the coordinates of the 1st corner, to plot a nice closed loop

        x = append(cornersXmm, cornersXmm[0])
        y = append(cornersYmm, cornersYmm[0])
       
        ax.plot(x, y, c=color[ccdCode])

        # Overplot the row closest to the readout register with a thicker line

        ax.plot([x[0], x[1]], [y[0], y[1]], c=color[ccdCode], linewidth=4)

        # If required, also plot the CCD labels

        if plotCCDlabels:

            minX = np.min(cornersXmm)
            maxX = np.max(cornersXmm)
            minY = np.min(cornersYmm)
            maxY = np.max(cornersYmm)
            middleX = minX + (maxX - minX) / 2.
            middleY = minY + (maxY - minY) / 2.
            ax.text(middleX, middleY, ccdCode, fontsize=45, color="gray")


    ax.set_xlabel("xFP [mm]")
    ax.set_ylabel("yFP [mm]")

    # That's it

    return










def drawSubfieldInFocalPlane(ccdCode, xCCD, yCCD, subfieldSizeX, subfieldSizeY, pixelSize):

    """
    PURPOSE: Draw a subfield in the focal plane.
    
             
    INPUT:   xCCD:          center x coordinate of the subfield [pixels]
             yCCD:          center y coordinate of the subfield [pixels]
             subfieldSizeX: size of the subfield along the x-axis [pixels]
             subfieldSizeY: size of the subfield along the y-axis [pixels]
             pixelSize:     the size of a pixel in microns
    
    OUTPUT:  Subfield frames draw on the current plot. The subfield is drawn with respect of the 
             coordinate system of the CCD.
    
             A blue dot indicates the lower left corner of the subfield.
             A green dot indicates the upper right corner of the subfield.
             A red dot indicates the center of the subfield.
    """

    # Compute the position of the subfield in pixel coordinates, for the current CCD, 
    # disregarding the physical extend of the CCD. LL = lower left, UR = upper right.

    zeroPointXmm = CCD[ccdCode]["zeroPointXmm"]
    zeroPointYmm = CCD[ccdCode]["zeroPointYmm"]
    ccdAngle     = CCD[ccdCode]["angle"]

    xFPprime, yFPprime = pixelToFocalPlaneCoordinates(xCCD, yCCD, pixelSize, zeroPointXmm, zeroPointYmm, ccdAngle)
    xFPprimeLL, yFPprimeLL = pixelToFocalPlaneCoordinates(xCCD - subfieldSizeX/2, yCCD - subfieldSizeY/2, \
                                                          pixelSize, zeroPointXmm, zeroPointYmm, ccdAngle)
    xFPprimeUR, yFPprimeUR = pixelToFocalPlaneCoordinates(xCCD + subfieldSizeX/2, yCCD + subfieldSizeY/2, \
                                                          pixelSize, zeroPointXmm, zeroPointYmm, ccdAngle)


    verts = [
        (xFPprimeLL, yFPprimeLL), # left, bottom
        (xFPprimeLL, yFPprimeUR), # left, top
        (xFPprimeUR, yFPprimeUR), # right, top
        (xFPprimeUR, yFPprimeLL), # right, bottom
        (xFPprimeLL, yFPprimeLL), # ignored
    ]

    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]

    path = Path(verts, codes)

    # Get the current axis

    currentAxis = plt.gca()
    currentAxis.add_patch( patches.PathPatch(path, facecolor='none') )

    plt.plot(xFPprime, yFPprime, 'ro')
    plt.plot(xFPprimeLL, yFPprimeLL, 'bo')
    plt.plot(xFPprimeUR, yFPprimeUR, 'go')
    plt.draw()


    return










def drawStarInFocalPlane(sim, raStar, decStar):

    """
    PURPOSE:  Draw a star given by the equatorial coordinates in the focal plane.

    INPUT:    sim:     instance of simulation class (see simulation.py)
              raStar:  right ascension of the star [rad]
              decStar: declination of the star     [rad]
    
    OUTPUT:   Draw a red dot where the star is located on the CCD

    TODO: Update doc-string

    """

    normal = True  # FIXME: where can we specify that we use the fast or normal Camera


    if sim["PSF/Model"] == "MappedFromFile":
        includeFieldDistortion = True
        isMapped               = True
        pathToPsfFile          = sim["PSF/MappedFromFile/Filename"]
        distortionCoefficients = None
    if (sim["Camera/IncludeFieldDistortion"] == "yes")  or (sim["Camera/IncludeFieldDistortion"] == "1"):
        includeFieldDistortion = True
        isMapped               = False
        pathToPsfFile          = None 
        FIELD_DISTORTION["Coeff"] = sim["Camera/FieldDistortion/Coefficients"]
        FIELD_DISTORTION["InverseCoeff"] = sim["Camera/FieldDistortion/InverseCoefficients"]
        distortionCoefficients = sim["Camera/FieldDistortion/Coefficients"]
    else:
        includeFieldDistortion = False
        distortionCoefficients = None
        pathToPsfFile          = None
        isMapped               = False 
        

    pixelSize             = float(sim["CCD/PixelSize"])
    raPlatform            = np.radians(float(sim["ObservingParameters/RApointing"]))
    decPlatform           = np.radians(float(sim["ObservingParameters/DecPointing"]))
    solarPanelOrientation = np.deg2rad(float(sim["Platform/SolarPanelOrientation"]))
    azimuthTelescope      = np.deg2rad(float(sim["Telescope/AzimuthAngle"]))
    tiltTelescope         = np.deg2rad(float(sim["Telescope/TiltAngle"]))
    focalPlaneAngle       = np.radians(float(sim["Camera/FocalPlaneOrientation/ConstantValue"]))
    focalLength           = float(sim["Camera/FocalLength/ConstantValue"]) * 1000.0  # [m] -> [mm]
    ccdZeroPointX         = float(sim["CCD/OriginOffsetX"])
    ccdZeroPointY         = float(sim["CCD/OriginOffsetY"])
    ccdAngle              = np.radians(float(sim["CCD/Orientation"]))

    xFPmm, yFPmm = skyToFocalPlaneCoordinates(raStar, decStar, raPlatform, decPlatform, solarPanelOrientation, tiltTelescope, azimuthTelescope, \
                                              focalPlaneAngle, focalLength)

    if includeFieldDistortion:
        if isMapped:
            xFPmm, yFPmm = mappedUndistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm, pathToPsfFile)
        else:
            xFPmm, yFPmm = undistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm, distortionCoefficients, focalLength)

    #ccdCode, xCCD, yCCD = getCCDandPixelCoordinates(raStar, decStar, raPlatform, decPlatform, tiltTelescope, azimuthTelescope,  \
    #                                                focalPlaneAngle, focalLength, pixelSize, includeFieldDistortion, FIELD_DISTORTION["Coeff"], normal)
    ccdCode, xCCD, yCCD = getCCDandPixelCoordinates(raStar, decStar, raPlatform, decPlatform, solarPanelOrientation, tiltTelescope, azimuthTelescope, focalPlaneAngle, focalLength, pixelSize, includeFieldDistortion, normal, isMapped,  distortionCoefficients, pathToPsfFile)

    if ccdCode == None:
        print ("Warning: DrawStarInFocalPlane(): The star doesn't fall on any of the CCDs.")
    else:
        drawPixelInFocalPlane(ccdCode, xCCD, yCCD, pixelSize)

    return













def drawPixelInFocalPlane(ccdCode, xCCD, yCCD, pixelSize):

    """
    PURPOSE: Plot a pixel from a particular CCD in the focal plane. The actual position in millimeter
             is shown as a red dot, while the pixel itself is drawn as a rectangle with edge pixelSize.

    INPUTS:  ccdCode:   for nominal camera: either '1', '2', '3', or '4'
                        for fast camer: either '1F', '2F', '3F', '4F'
             xCCDpix:   x-coordinate (column number, zero-based) of the pixel on the CCD  [pix]
             yCCDpix:   y-coordinate (row number, zero-based) of the pixel on the CCD  [pix]
             pixelSize: the size of a pixel in micron

    OUTPUTS: None
    """

    # Compute the position of the star in pixel coordinates, for the current CCD, 
    # disregarding the physical extend of the CCD

    zeroPointXmm = CCD[ccdCode]["zeroPointXmm"]
    zeroPointYmm = CCD[ccdCode]["zeroPointYmm"]
    ccdAngle     = CCD[ccdCode]["angle"]

    xFPmm, yFPmm = pixelToFocalPlaneCoordinates(xCCD, yCCD, pixelSize, zeroPointXmm, zeroPointYmm, ccdAngle)   # [mm]

    # Get the current axis

    currentAxis = plt.gca()
    currentAxis.add_patch( patches.Rectangle( (xFPmm, yFPmm), pixelSize / 1000.0, pixelSize / 1000.0, fill=False) )

    plt.plot(xFPmm, yFPmm, 'ro')
    plt.draw()
    plt.show()

    return








def skyProjection(longitude, latitude, fig, origin=0, projection="mollweide"):
    
    """
    Plot sources with coordinates longitude & latitude (equatorial, galactic,...)
    in the projected sky.
    
    @param longitude:  longitude coordinate in [0, 360]                      [deg]
    @param latitude:   latitude coordinate in [-90, +90]                     [deg]
    @param fig:        matplotlib figure, the output of plt.figure()
    @param origin:     longitude value in the center of the plot             [deg]
    @param projection: either 'mollweide', 'aitoff', 'hammer', or 'lambert'
    
    @return axes: the output of plt.scatter() in the given figure.
    
    """
    
    # Shift the longitude values around the origin
    
    longitude = np.remainder(longitude+360-origin, 360) 
    
    # Rescale the range from [0,360] to [-180, +180]
    
    longitude[longitude>180] -=360

    # Reverse the longitude so that East is to the left
    
    longitude = -longitude
    
    # Create the plot
    
    axes = fig.add_subplot(111, projection=projection)
    axes.grid(True)
    
    # Adapt the tick labels on the x-axis to take into account the origin shift
    
    tickLabels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tickLabels = np.remainder(tickLabels+360+origin, 360)
    axes.set_xticklabels(tickLabels) 
 
    # Plot the sources
    
    axes.scatter(np.radians(longitude), np.radians(latitude), s=3)
    
    # That's it!
    
    return axes







def drawStarsInSkyAitoff(fig, raStars, decStars, magStars, skymap=None, cbarOrientation=None, cbarMap='rainbow'):
    """
    Project and plot a catalog of stars on the sky in a Aitoff Galactic projection.
    This plot uses the astropy library to make the ICRS to Galactic coordinate
    transformation together with a nice Galactic background image. To show the plot
    it is necessary to introduce a "plt.show()" after the function call. This module
    scales the scatter plot markersize of the stars according to their sample size.

    Parameters
    ----------
    fig : object
        Figure matplotlib.pyplot object to define e.g. figsize
    raStars : list, array
        Right ascension of stars [deg]
    decStars : list, array
        Declination of stars [deg]
    magStars : list, array
        Magnitudes of stars
    cbarOrientation : str
        Colorbar orientation. Default 'horizontal' else 'vertical'
    cbarMap : str
        Colormap of colorbar. Default 'rainbow'

    Return
    ------
    axes : object
        Axes matplotlib.pyplot handle object to be modified by the user
    """

    # Convert coordinates from ICRS to Galactic using astropy

    gal = SkyCoord(raStars, decStars, frame='icrs', unit=u.deg)
    gal = gal.galactic

    # Plot Aitoff projection in Galactic coordinates

    plt.title('Aitoff projection in Galactic coordinates', fontsize=18, y=1.02)
    fig, ax = fig
    fs = 16
    if len(raStars) <= 1e2: ms = 3.
    if len(raStars) >= 1e2 and len(raStars) < 1e3: ms = 1.3
    if len(raStars) >= 1e3 and len(raStars) < 1e4: ms = 1.
    if len(raStars) >= 1e4: ms = 0.1

    # Plot Galactic map as background (e.g. Gaia DR3)
    # E.g.: skymap = plt.imread('skymap.png')

    if skymap is not None:
        ax.imshow(skymap)

    # Add the sky projection ontop as transparent layer

    axes = fig.add_subplot(111, projection='aitoff', facecolor='none')

    # Plot the targets on the sky (autumn_r, rainbow)

    im = plt.scatter(-gal.l.wrap_at('180d').radian, gal.b.radian, c=magStars, s=ms, cmap=cbarMap, zorder=3)

    # Vertical or horizontal colorbar showing magnitudes

    if cbarOrientation == 'vertical':
        cbarax = fig.add_axes([0.905, 0.2, 0.02, 0.57])
        cbar = plt.colorbar(im, orientation='vertical', cax=cbarax, extend='both')
        cbar.set_label(r'PLATO passband, $P$', fontsize=fs)
        cbar.ax.tick_params(labelsize=fs)
    else:
        cbarax = fig.add_axes([0.25, 0.06, 0.525, 0.03])
        cbar = plt.colorbar(im, orientation='horizontal', cax=cbarax, extend='both')
        cbar.set_label(r'PLATO passband, $P$', fontsize=fs)
        cbar.ax.tick_params(labelsize=fs)

    # Change the tick labels so that they are 0->360, rather than -180->+180

    tickLabels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tickLabels = np.remainder(tickLabels+360, 360)
    axes.set_xticklabels(tickLabels)

    # Change y ticks and remove last to make space for title

    tickLabels = np.array([-75, -60, -45, -30, -15, 0, 15, 30, 45, 60, ''])
    axes.set_yticklabels(tickLabels)

    # Change color of x tick labels

    axes.tick_params(axis='x', colors='w')

    # Increase x and y tick labels

    axes.xaxis.set_tick_params(labelsize=fs+1)
    axes.yaxis.set_tick_params(labelsize=fs)

    # Add axis labels

    axes.set_xlabel(r'Longitude, $l$ [deg]', fontsize=fs)
    axes.set_ylabel(r'Latitude, $b$ [deg]', fontsize=fs)

    # Set grid and remore outer ticks (if set by default)

    axes.grid(True, alpha=0.3)
    ax.axis('off')
    plt.draw()

    # That's it

    return axes









def plotStellarSampleDistributions(fig, mag, magCon, magRange, numConPerTar, distCon):
    """
    This function plots 4 different stellar sample distribution plots
    for an PLATO Input Catalogue (PIC)
    1) Magnitude distribution of PIC targets
    2) Magnitude distribution of PIC contaminants
    3) Number distribution of contaminants per target
    4) Distance distribution of contaminants

    Parameters
    ----------
    mag : list, array
        The stellar target magnitudes
    magCon : list, array
        The stellar contaminant magnitudes
    magRange : list, array
        Upper and lower magnitude limit for input sample
    numConPerTar : list, array
        The number of contaminants (integer) per target star
    distCon : list, array
        Distances of each contaminant star w.r.t. their target

    Return
    ------
    axes : object
        Axes matplotlib.pyplot handle object to be modified by the user
    """

    # Copy figure object not to overwrite it

    fig1, axes = fig

    # Prepare bins and plot magnitude distribution of targets

    magbinTar  = 0.1
    if magRange[1]-magRange[0] > 10: magbinTar = 0.2
    binsizeTar = int((magRange[1] - magRange[0]) / magbinTar) + 1
    binlistTar = np.linspace(magRange[0], magRange[1], binsizeTar)

    axes[0,0].hist(mag, binlistTar, facecolor='b', edgecolor='b', fill=True, alpha=0.3)
    axes[0,0].set_title('Magnitude distribution of PIC targets')
    axes[0,0].set_xlabel('Gaia V Magnitude')
    axes[0,0].set_ylabel('Number of stars')
    axes[0,0].locator_params(axis='y', integer=True)
    axes[0,0].tick_params(axis='x', which='minor', bottom=True, top=False)
    axes[0,0].tick_params(axis='x', which='major', bottom=True, top=False)
    axes[0,0].tick_params(axis='y', which='minor', left=False, right=False)
    axes[0,0].tick_params(axis='y', which='major', left=True, right=False)
    axes[0,0].grid(axis='y', color='gray', alpha=0.3)

    # Prepare bins and plot magnitude distribution of contaminants

    magbinCon  = 0.2
    binsizeCon = int((np.max(magCon) - np.min(magCon)) / magbinCon) + 1
    binlistCon = np.linspace(round(np.min(magCon)), round(np.max(magCon)), binsizeCon)

    axes[0,1].hist(magCon, binlistCon, facecolor='m', edgecolor='m', fill=True, alpha=0.3)
    axes[0,1].set_title('Magnitude distribution of PIC contaminants')
    axes[0,1].set_xlabel('Gaia V Magnitude')
    axes[0,1].set_ylabel('Number of stars')
    axes[0,1].tick_params(axis='x', which='minor', bottom=True, top=False)
    axes[0,1].tick_params(axis='x', which='major', bottom=True, top=False)
    axes[0,1].tick_params(axis='y', which='minor', left=False, right=False)
    axes[0,1].tick_params(axis='y', which='major', left=True, right=False)
    axes[0,1].grid(axis='y', color='gray', alpha=0.3)

    # Prepare bins and plot number distribution of contaminants per target

    numbinCon  = 1
    binsizeNum = int((np.max(numConPerTar) - 0) / numbinCon) + 2
    binlistNum = np.linspace(-0.5, np.max(numConPerTar)+0.5, binsizeNum)  # -0.5 because num x-axis

    axes[1,0].hist(numConPerTar, binlistNum, facecolor='g', edgecolor='g', fill=True, log=True, alpha=0.3)
    axes[1,0].set_title('Number distribution of contaminants per target')
    axes[1,0].set_xlabel('Number of contaminants')
    axes[1,0].set_ylabel('Number of targets')
    axes[1,0].tick_params(axis='x', which='minor', bottom=False, top=False)
    axes[1,0].tick_params(axis='x', which='major', bottom=True, top=False)
    axes[1,0].tick_params(axis='y', which='minor', left=False, right=False)
    axes[1,0].tick_params(axis='y', which='major', left=True, right=False)
    axes[1,0].grid(axis='y', color='gray', alpha=0.3)

    # Prepare bins and plot distance distribution of contaminants in respect to their target star

    distbinCon  = 1.0
    binsizeDist = int((np.max(distCon) - np.min(distCon)) / distbinCon) + 2  # +1 extra because zero is rare
    binlistDist = np.linspace(round(np.min(distCon)), round(np.max(distCon)), binsizeDist)

    axes[1,1].hist(distCon, binlistDist, facecolor='orange', edgecolor='orange', fill=True, alpha=0.4)
    axes[1,1].set_title('Distance distribution of contaminants')
    axes[1,1].set_xlabel('Distances [arcsec]')
    axes[1,1].set_ylabel('Number of stars')
    axes[1,1].locator_params(axis='y', integer=True)
    axes[1,1].tick_params(axis='x', which='minor', bottom=True, top=False)
    axes[1,1].tick_params(axis='x', which='major', bottom=True, top=False)
    axes[1,1].tick_params(axis='y', which='minor', left=False, right=False)
    axes[1,1].tick_params(axis='y', which='major', left=True, right=False)
    axes[1,1].grid(axis='y', color='gray', alpha=0.3)

    # Layout

    plt.tight_layout()

    # That's it!

    return axes









def plotYawPitchRollTimeSeries(fig, time, signal, units, title=False, ylims=False):
    """
    Function to plot the time series of yaw, pitch, and roll for both AOSC jitter and thermo drift.
    Along with the time series plots the Root-Mean-Square (RMS) are calculated and plotted in each
    planel, respectively.

    Parameters
    ----------
    time : ndarray
        Array of time points. Make sure time units matches label.
    signals : ndarray, list-ndarray
        List or array of Yaw, Pitch, and Roll time series. Make sure that ampltide unit matches label.
    units : list-str
        List of strings of physical unit: ['time-unit', 'ampltude-unit']
    title : str (optional)
        Title in a string
    ylims : list-str (optional)
        List with ymin and ymax limits

    Return
    ------
    axes : object
        Axes matplotlib.pyplot handle object to be modified by the user
    """

    # Datasets to loop over

    numData = len(signal)

    # Handle yaxis limits

    if ylims is False:
        lim = 0.5*np.max(np.abs(signal))

    # Adjust linewidth after data

    if len(time) < 1e3: lw = 1
    else: lw = 0.5

    # Make plot

    labels = ['Yaw', 'Pitch', 'Roll']
    colors = ['royalblue', 'lightseagreen', 'limegreen']

    #fig, axes = plt.subplots(numData, 1, figsize=figsize)

    for plot in range(numData):

        axes = fig.add_subplot(numData, 1, plot+1)

        # Make sure that time series is redual around zero

        signal[plot] -= np.median(signal[plot])

        # Plot timeseries

        axes.plot(time, signal[plot], '-', c=colors[plot], lw=lw)

        # Add root-mean-square lines

        rms = np.sqrt(np.mean(signal[plot]**2))
        axes.axhline(+rms, c='k', ls='--', lw=0.7, label='RMS = {0:.3f} {1}'.format(rms, units[1]))
        axes.axhline(-rms, c='k', ls='--', lw=0.7)
        axes.legend(loc='upper right')

        # Latter settings

        axes.set_ylabel('{0} [{1}]'.format(labels[plot], units[1]))
        axes.set_xlim(np.min(time), np.max(time))
        axes.set_ylim(-lim, +lim)

        # Remove tick labels on x axis except for last plot

        if plot < numData-1:
            axes.tick_params(labelbottom=False)
            #axes.set_ylim(-lim, +lim)
        else:
            axes.set_xlabel('Time [{0}]'.format(units[0]))

        # Title

        if plot == 0: axes.set_title(title, fontsize=fs)

    # Remaining

    #if title is not False: fig.text(0.5, 0.95, title, ha='center', fontsize=fs)
    plt.tight_layout()
    plt.subplots_adjust(hspace = .001)

    # Finito!

    return axes









def plotYawPitchRollPSD(fig, time, signals, xmin=False, ylim=False, carbox=False, title=False, labels=False, misreq=False):
    """
    This function takes a Yaw, Pitch, and Roll time series and plots the Power Spectral Density (PSD)
    function for each. Alongside the data a median filter is plotted with a default carbox length of
    144 time points, corresponding to 1 hour precision if the time series are given in seconds.

    Parameters
    ----------
    time : narray
        Time points [s]
    signals : narray, list-narray
        Either single signal array or a list of signal arrays
    xmin : float (optional)
        Limit for x min. The x max limit is the Nyquist frequency
    ylim : list-float (optional)
        List of y min and max limit ["y-min", "y-max"]
    carbox : int (optional)
        Length of median carbox filter. Default is 3600s/25s = 144
    title : str (optional)
        Title for plot
    labels : list-str (optinal)
        List of string labels where the first is the xlabel and the rest is ylabels
    misreq : bool (optional)
        If "True" the mission requirements for the AOSC will be plotted alonside the data.

    Return
    ------
    Plot or/and saved plot to PNG.
    """

    # Number of data sets

    numData = len(signals)

    # Find time step

    sampling = np.diff(time)[0]

    # Choose carbox length of 1 hour if time is in seconds

    if carbox is False: carbox = 144

    # Make plot

    labels = ['Yaw', 'Pitch', 'Roll']
    colors = ['tomato', 'darkorange', 'gold']

    for plot in range(numData):

        # Create axes objects

        axes = fig.add_subplot(numData, 1, plot+1)

        # Find PSD and median filter

        freq, PSD = powerDensityFFT(signals[plot], sampling)
        PSD_med   = median_filter(PSD, carbox)
        #freq      *= 1e-3

        # Plot results

        axes.plot(freq, PSD,     '-', c=colors[plot], lw=lw, label=labels[plot])
        axes.plot(freq, PSD_med, 'k-', lw=lw+1)

        # Plot mission requirements

        if misreq:
            axes.plot([1e-1, 2e1], [1e10, 1e1], c='k', linestyle='--', lw=1)
            axes.plot([2e1,  1e4], [1e1,  1e1], c='k', linestyle='--', lw=1)

        # Log scaling

        axes.set_xscale("log")
        axes.set_yscale("log")

        # Latter settings

        axes.set_ylabel('{0}'.format(labels[plot]) + r' [arcsec$^2$ Hz$^{-1}$]')

        # Remove tick labels on x axis except for last plot

        if plot < numData-1:
            axes.tick_params(labelbottom=False)

        # Set x-min limit

        if xmin is not False: axes.set_xlim(xmin, freq.max())
        #else: axes.set_xlim(freq.min(), freq.max())

        # Set y limits

        if ylim is not False:
            axes.set_ylim(ylim[0], ylim[1])


        # Remove tick labels on x axis except for last plot

        if plot < numData-1: axes.tick_params(labelbottom=False)

        # Set legends

        axes.legend(loc='upper right')

        # Set title

        if title is not False and plot == 0: axes.set_title(title, fontsize=fs)

    # Remaining

    plt.xlabel(r'Frequency [mHz]')
    plt.tight_layout()
    plt.subplots_adjust(hspace = .001)

    # Finito!

    return axes










def plotYawPitchRollJitter(time, signals, clabel, cmap='gnuplot', plottype='short', tpoint=100, title=False):
    """
    This function can be used to plot a time series of spacecraft jitter.
    For time series on short time scales, the correlation between yaw, pitch,
    and roll can be illustrated using the the "plottype='short'" option. For
    visualising the entire jitter time series correlated between yaw, pitch,
    and rool use instead the "plottype='long'" option.

    Parameters
    ----------
    time : ndarray
        Array of times.
    signals : list, ndarray
        List or array of individual signal arrays, i.e. [yaw, pitch, roll]
    clabel : str
        String with the color-bar time-label.
    cmap : str
        String specifying the matplotlib color map. Default is "gnuplot"
    plottype : str
        String determine the time scale of the plot. Options are "short" and "long"
    tpoint : int
        Number time points to use for the short time scale visualition only
    title : str
        String of the main title.

    Return
    ------
    Nothing so far
    """

    # Hardcode values

    lim    = 0.125  #np.max(np.abs(data))
    nticks = 5
    labels = ['Yaw [arcsec]', 'Pitch [arcsec]', 'Roll [arcsec]']
    sms = 5
    al  = 1

    # Adjust parameters to a zero-point level

    time = time - np.min(time)
    signals[0] -= np.median(signals[0])
    signals[1] -= np.median(signals[1])
    signals[2] -= np.median(signals[2])

    # PLOT CORRELATIONS FOR SHORT TERM JITTER

    if plottype == 'short':

        fig, ax = plt.subplots(3, 3, figsize=(10, 7))

        for row in range(3):

            # Plots

            ax[row, 0].plot(signals[1][tpoint*row:tpoint*(row+1)], signals[0][tpoint*row:tpoint*(row+1)], 'k-', alpha=al, lw=lw, zorder=1)
            ax[row, 1].plot(signals[2][tpoint*row:tpoint*(row+1)], signals[0][tpoint*row:tpoint*(row+1)], 'k-', alpha=al, lw=lw, zorder=1)
            ax[row, 2].plot(signals[2][tpoint*row:tpoint*(row+1)], signals[1][tpoint*row:tpoint*(row+1)], 'k-', alpha=al, lw=lw, zorder=1)
            im0 = ax[row, 0].scatter(signals[1][tpoint*row:tpoint*(row+1)], signals[0][tpoint*row:tpoint*(row+1)], c=time[tpoint*row:tpoint*(row+1)], s=sms, cmap=cmap, zorder=2)
            im1 = ax[row, 1].scatter(signals[2][tpoint*row:tpoint*(row+1)], signals[0][tpoint*row:tpoint*(row+1)], c=time[tpoint*row:tpoint*(row+1)], s=sms, cmap=cmap, zorder=2)
            im2 = ax[row, 2].scatter(signals[2][tpoint*row:tpoint*(row+1)], signals[1][tpoint*row:tpoint*(row+1)], c=time[tpoint*row:tpoint*(row+1)], s=sms, cmap=cmap, zorder=2)

            # Labels

            ax[2, 0].set_xlabel(labels[1])
            ax[2, 1].set_xlabel(labels[2])
            ax[2, 2].set_xlabel(labels[2])
            ax[row, 0].set_ylabel(labels[0])
            ax[row, 1].set_ylabel(labels[0])
            ax[row, 2].set_ylabel(labels[1])

            # Remove tick labels on x axis except for last plot

            if row < 2:
                ax[row, 0].tick_params(labelbottom=False)
                ax[row, 1].tick_params(labelbottom=False)
                ax[row, 2].tick_params(labelbottom=False)

            # Duplicate settings for each plotted row

            for col, im in zip(range(3), [im0, im1, im2]):

                # Axes limits

                ax[row, col].set_xlim(-lim, lim)
                ax[row, col].set_ylim(-lim, lim)
                ax[row, col].set_aspect('equal', 'box')

                # Force the same number of ticks

                ax[row, col].xaxis.set_major_locator(MaxNLocator(nticks))
                ax[row, col].yaxis.set_major_locator(MaxNLocator(nticks))

                # Set grid

                ax[row, col].grid(c='gray', ls='-', lw=lw, alpha=al)

                # Color bars

                div = make_axes_locatable(ax[row, col])
                cax = div.append_axes('right', size='10%', pad=0.1)
                cbar = plt.colorbar(im, ax=ax[row, col], cax=cax, extend='max')
                cbar.ax.invert_yaxis()
                if im == im2:
                    cbar.set_label(clabel)
                else:
                    cbar.remove()


    # PLOT CORRELATIONS FOR ENTIRE TIMESERIES

    if plottype == 'long':

        # Limits and grid
        lim = 0.28  #np.max(np.abs(signals))
        nticks = 6
        time = time/(60*60)

        fig1, ax1 = plt.subplots(1, 3, figsize=(10, 2.8))

        # Plot
        ax1[0].plot(signals[1], signals[0], 'k-', alpha=al, lw=lw, zorder=1)
        ax1[1].plot(signals[2], signals[0], 'k-', alpha=al, lw=lw, zorder=1)
        ax1[2].plot(signals[2], signals[1], 'k-', alpha=al, lw=lw, zorder=1)
        im0 = ax1[0].scatter(signals[1], signals[0], c=time, s=2, cmap='magma', zorder=2)
        im1 = ax1[1].scatter(signals[2], signals[0], c=time, s=2, cmap='magma', zorder=2)
        im2 = ax1[2].scatter(signals[2], signals[1], c=time, s=2, cmap='magma', zorder=2)

        # Labels
        ax1[0].set_xlabel(labels[1])
        ax1[0].set_ylabel(labels[0])
        ax1[1].set_xlabel(labels[2])
        ax1[1].set_ylabel(labels[0])
        ax1[2].set_xlabel(labels[2])
        ax1[2].set_ylabel(labels[1])

        # Duplicate settings
        for plot in range(3):

            # Adjust axes
            ax1[plot].set_xlim(-lim, lim)
            ax1[plot].set_ylim(-lim, lim)
            ax1[plot].set_aspect('equal', 'box')

            # Force the same number of ticks
            ax1[plot].xaxis.set_major_locator(MaxNLocator(nticks))
            ax1[plot].yaxis.set_major_locator(MaxNLocator(nticks))

            # Plot grid
            ax1[plot].set_axisbelow(False)
            ax1[plot].grid(c='gray', ls='-', lw=lw, alpha=al)

            # We put the colorbar for then remove
            # This is need to keep the same size of each subplot..
            div1 = make_axes_locatable(ax1[2])
            cax1 = div1.append_axes("right", size="10%", pad=0.1)
            cbar1 = plt.colorbar(im2, ax=ax1[2], cax=cax1, extend='max')
            cbar1.ax.invert_yaxis()
            if plot == 2:
                cbar1.set_label('Time [hours]')
            else:
                cbar1.remove()

    # Adjust subplot spacing

    fig.tight_layout()
    fig.subplots_adjust(hspace = .001)
    plt.show()
    # Finito!

    return
