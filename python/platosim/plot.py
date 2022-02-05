#!/usr/bin/env python3

import os
import h5py
import numpy as np

from scipy import constants as c
from scipy.ndimage import median_filter

from matplotlib import pyplot as plt
from matplotlib.pyplot import cm
from matplotlib import patches
from matplotlib.path import Path
from matplotlib.ticker import MaxNLocator, ScalarFormatter
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib.animation as animation

import astropy.units as u
from astropy.coordinates import SkyCoord

from platosim.simfile import SimFile
import platosim.referenceFrames as rf
import platosim.utilities as ut
import platosim.noise as ns

# Top level Matplotlib settings to ease the writing

fs = 15    # Font size
lw = 0.3   # Line width
ms = 0.5   # Marker size

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








def drawCCDsInFocalPlane(pixelSize=18, plotCCDlabels=True, normal=True):

    """
    PURPOSE: Plot the 4 CCDs in the focal plane in the FP' reference frame.
             May serve as a background to overplot the projected stars on the focal plane

    INPUT: pixelSize: size of 1 pixel [micron]. Default 18 micron for PLATO.
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

        cornersXmm, cornersYmm = rf.computeCCDcornersInFocalPlane(ccdCode, pixelSize)

        # Repeat the coordinates of the 1st corner, to plot a nice closed loop

        x = np.append(cornersXmm, cornersXmm[0])
        y = np.append(cornersYmm, cornersYmm[0])

        ax.plot(x, y, c=color[ccdCode])

        # Overplot the row closest to the readout register with a thicker line

        ax.plot([x[0], x[1]], [y[0], y[1]], c=color[ccdCode], linewidth=4)

        # If required, also plot the CCD labels

        if plotCCDlabels:

            minX = np.min(cornersXmm)
            maxX = np.max(cornersXmm)
            minY = np.min(cornersYmm)
            maxY = np.max(cornersYmm)
            middleX = minX + (maxX - minX) / 8.
            middleY = minY + (maxY - minY) / 8.
            ax.text(middleX, middleY, ccdCode, fontsize=30, color="gray")

    ax.set_xlabel('$x_{FP}$ [mm]', fontsize = 20)
    ax.set_ylabel('$y_{FP}$ [mm]', fontsize = 20)

    # Finito!

    return











def drawCCDsInCameraFocalPlane(fig):
    """
    PURPOSE: Draw the CCDs in the focal plane of a N-CAM.

    INPUT:   fig:  A matplotlib.pyplot figure object [e.g. plt.figure(figsize=(10,10)].

    OUTPUT:  ax:   Axes object to modify the figure object.
    """

    # Constants

    pixelSize   = 18        # [µm]
    fovDegrees  = 18.8908   # [deg]
    focalLength = 247.52    # [mm]
    ccdCodes    = ["1", "2", "3", "4"]

    # Size of the FOV

    fovMm = focalLength * np.tan(np.radians(fovDegrees))

    # Star plotting

    ax = fig.add_subplot(111)

    # Plot the camera aperture

    circ = plt.Circle((0, 0), radius=fovMm, color="none", linewidth=2, label="Telescope FOV", zorder=1)
    ax.add_patch(circ)
    circ.set_edgecolor("g")
    circ.set_facecolor("lightgray")

    # Plot location of optical axis

    ax.plot([0], [0], "rx", label="Optical axis", zorder=2)

    # Plot CCDs

    ori = [[-30, 0, 0, -30], [30, 0, 0, -30], [30, 0, 0, 30], [-30, 0, 0, 30]]

    for ccdCode in ccdCodes:

        # Fetch CCD corners

        cornersX, cornersY = rf.computeCCDcornersInFocalPlane(ccdCode, pixelSize)

        # Draw CCD names at the center of each CCD

        ax.text(np.mean(cornersX) - 10, np.mean(cornersY), "CCD " + ccdCode, fontsize=20)

        # Draw each CCD

        cornersX = np.append(cornersX, cornersX[0])
        cornersY = np.append(cornersY, cornersY[0])

        if ccdCode == "1":
            ax.plot(cornersX, cornersY, color="b", label="CCD footprint", zorder=3)
        else:
            ax.plot(cornersX, cornersY, color="b", zorder=3)

        # Plot arrays to indicate CCD origin

        ax.arrow(cornersX[0], cornersY[0], ori[int(ccdCode)-1][0], ori[int(ccdCode)-1][1],
                 head_width=3, head_length=3, fc='k', ec='k', linewidth=2, zorder=4)
        ax.arrow(cornersX[0], cornersY[0], ori[int(ccdCode)-1][2], ori[int(ccdCode)-1][3],
                 head_width=3, head_length=3, fc='k', ec='k', linewidth=2, zorder=4)

    # Settings

    plt.legend(prop={'size': 20}, bbox_to_anchor=(1.0, 1.0))
    ax.set_title('CCDs in camera focal plane', fontsize = 24)
    ax.set_xlabel(r'$x_{\mathrm{FP}}$ [mm]', fontsize = 20)
    ax.set_ylabel(r'$y_{\mathrm{FP}}$ [mm]', fontsize = 20)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)

    # Finito!

    return ax








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

    zeroPointXmm = rf.CCD[ccdCode]["zeroPointXmm"]
    zeroPointYmm = rf.CCD[ccdCode]["zeroPointYmm"]
    ccdAngle     = rf.CCD[ccdCode]["angle"]

    xFPprime, yFPprime = rf.pixelToFocalPlaneCoordinates(xCCD, yCCD, pixelSize, zeroPointXmm, zeroPointYmm, ccdAngle)
    xFPprimeLL, yFPprimeLL = rf.pixelToFocalPlaneCoordinates(xCCD - subfieldSizeX/2, yCCD - subfieldSizeY/2, \
                                                             pixelSize, zeroPointXmm, zeroPointYmm, ccdAngle)
    xFPprimeUR, yFPprimeUR = rf.pixelToFocalPlaneCoordinates(xCCD + subfieldSizeX/2, yCCD + subfieldSizeY/2, \
                                                             pixelSize, zeroPointXmm, zeroPointYmm, ccdAngle)


    verts = [
        (xFPprimeLL, yFPprimeLL),  # left, bottom
        (xFPprimeLL, yFPprimeUR),  # left, top
        (xFPprimeUR, yFPprimeUR),  # right, top
        (xFPprimeUR, yFPprimeLL),  # right, bottom
        (xFPprimeLL, yFPprimeLL),  # ignored
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

    xFPmm, yFPmm = rf.skyToFocalPlaneCoordinates(raStar, decStar, raPlatform, decPlatform, solarPanelOrientation,
                                                 tiltTelescope, azimuthTelescope, focalPlaneAngle, focalLength)

    if includeFieldDistortion:
        if isMapped:
            xFPmm, yFPmm = rf.mappedUndistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm, pathToPsfFile)
        else:
            xFPmm, yFPmm = rf.undistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm, distortionCoefficients, focalLength)

    #ccdCode, xCCD, yCCD = getCCDandPixelCoordinates(raStar, decStar, raPlatform, decPlatform, tiltTelescope, azimuthTelescope,  \
    #                                                focalPlaneAngle, focalLength, pixelSize, includeFieldDistortion, FIELD_DISTORTION["Coeff"], normal)
    ccdCode, xCCD, yCCD = rf.getCCDandPixelCoordinates(raStar, decStar, raPlatform, decPlatform, solarPanelOrientation, tiltTelescope, azimuthTelescope, focalPlaneAngle, focalLength, pixelSize, includeFieldDistortion, normal, isMapped,  distortionCoefficients, pathToPsfFile)

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

    zeroPointXmm = rf.CCD[ccdCode]["zeroPointXmm"]
    zeroPointYmm = rf.CCD[ccdCode]["zeroPointYmm"]
    ccdAngle     = rf.CCD[ccdCode]["angle"]

    xFPmm, yFPmm = rf.pixelToFocalPlaneCoordinates(xCCD, yCCD, pixelSize, zeroPointXmm, zeroPointYmm, ccdAngle)   # [mm]

    # Get the current axis

    currentAxis = plt.gca()
    currentAxis.add_patch( patches.Rectangle( (xFPmm, yFPmm), pixelSize / 1000.0, pixelSize / 1000.0, fill=False) )

    plt.plot(xFPmm, yFPmm, 'ro')
    plt.draw()
    plt.show()

    return













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

    if normal:
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

        cornersXmm, cornersYmm = rf.computeCCDcornersInFocalPlane(ccdCode, pixelSize)

        # Compute the equatorial sky coordinates [rad] from the the focal plane FP' coordinates [mm] of the corners

        ra, dec = rf.focalPlaneToSkyCoordinates(cornersXmm, cornersYmm, raPlatform, decPlatform, solarPanelOrientation,
                                                tiltAngle, azimuthAngle, focalPlaneAngle, focalLength)

        # Repeat the coordinates of the 1st corner, to plot a nice closed loop
        # Convert from radians to degrees

        ra  = np.append(ra, ra[0])
        dec = np.append(dec, dec[0])

        # The sky projection expects a longitude in [-pi, +pi] rather than [0, 2* pi]
        # Moreover, the longitude should be reversed so that East is to the left

        ra[ra>np.pi] -= 2 * np.pi
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
        raRadians.append(-ra[index] * np.pi / 180.0)
        decRadians.append(dec[index] * np.pi / 180.0)

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














def skyProjection(fig, longitude, latitude, origin=0, projection="mollweide"):
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

    longitude[longitude > 180] -= 360

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











def plotPlatoFOV(pointingField, raStars, decStars, magStars=None, nCamVis=None, skymap=None, save=False):
    """

    Parameters
    ----------

    Return
    ------
    axes : object
        Axes matplotlib.pyplot handle object to be modified by the user
    """

    import ligo.skymap.plot

    # Select field
    
    indir = os.getenv('PLATO_WORKDIR') + '/platonium/pic/PIC1.1.0/'

    if pointingField == 'NPF': PF_gal = [65.0, 30.0]
    if pointingField == 'SPF': PF_gal = [253.0, -30.0]

    PF_gal  = SkyCoord(PF_gal[0], PF_gal[1], frame='galactic', unit='deg')  # [deg]
    PF_icrs = PF_gal.icrs  # [deg]

    PF06 = np.load(indir + f'{pointingField}-NCAM06.npy')
    PF12 = np.load(indir + f'{pointingField}-NCAM12.npy')
    PF18 = np.load(indir + f'{pointingField}-NCAM18.npy')
    PF24 = np.load(indir + f'{pointingField}-NCAM24.npy')

    starPF06 = SkyCoord(PF06[:,0]*u.deg, PF06[:,1]*u.deg, frame='icrs', unit='deg')
    starPF12 = SkyCoord(PF12[:,0]*u.deg, PF12[:,1]*u.deg, frame='icrs', unit='deg')
    starPF18 = SkyCoord(PF18[:,0]*u.deg, PF18[:,1]*u.deg, frame='icrs', unit='deg')
    starPF24 = SkyCoord(PF24[:,0]*u.deg, PF24[:,1]*u.deg, frame='icrs', unit='deg')
    
    # Load brightest stars
    starPF = SkyCoord(raStars*u.deg, decStars*u.deg, frame='icrs', unit='deg')

    # MAKE PLOTS
    fig = plt.figure(figsize=(10,10))
    ax = plt.axes(projection='astro zoom', center=PF_icrs, radius='35 deg', rotate='180 deg')

    # Plot PIC1.1.0 stars after N-CAM visibility

    ax.plot(starPF06.ra.deg, starPF06.dec.deg, '.', c='skyblue',
            transform=ax.get_transform('world'), markersize=1, zorder=1)
    ax.plot(starPF12.ra.deg, starPF12.dec.deg, '.', c='deepskyblue',
            transform=ax.get_transform('world'), markersize=1, zorder=2)
    ax.plot(starPF18.ra.deg, starPF18.dec.deg, '.', c='dodgerblue',
            transform=ax.get_transform('world'), markersize=1, zorder=3)
    ax.plot(starPF24.ra.deg, starPF24.dec.deg, '.', c='royalblue',
            transform=ax.get_transform('world'), markersize=1, zorder=4)

    # Plot saturated stars and add legend scaled to the stellar magnitudes

    if magStars is not None:
        maxMarkerSize = 30
        dm = (max(magStars) - magStars) * maxMarkerSize
        mag_range = np.arange(min(magStars), max(magStars)).astype(int)
        dm_range  = (max(magStars) - mag_range) * maxMarkerSize/10
        mark, color = 'o', 'gold'
        handle = [plt.plot([],[], "o", c='gray', ms=dm_range[i], ls="")[0] for i in range(len(dm_range))]
        ax.legend(handles=handle, labels=mag_range.tolist(), loc='upper right', title=r"P [mag]", fontsize=16, title_fontsize=16)
    else:
        dm, mark, color = 20, '*', 'none'
    # Plot all stars
    scatter = ax.scatter(starPF.ra.deg, starPF.dec.deg, transform=ax.get_transform('world'), 
                         s=dm, marker=mark, c=color, ec='k', lw=1, zorder=5)

    # Plot pointing of each camera group

    #camPointing = coorCameraGroup(PF_icrs.ra.deg, PF_icrs.dec.deg)
    #ax.plot(camPointing.ra.deg, camPointing.dec.deg, 'rx', transform=ax.get_transform('world'), markersize=10, zorder=6)

    # Plot pointing of PIC1.1.0 and PIC2.0.0

    ax.plot(PF_icrs.ra.deg, PF_icrs.dec.deg, '*', transform=ax.get_transform('world'), ms=20, c='k', mfc='r', zorder=6)
    ax.plot(277.18, 52.85, '*', transform=ax.get_transform('world'), ms=20, c='k', mfc='b', zorder=7)

    # Add-on's

    ax.scalebar((0.05, 0.05), 10 * u.deg).label()
    ax.compass(0.95, 0.05, 0.1)
    ax.grid(color='gray')

    # Settings

    ax.set_title(f'{pointingField} pointing', fontsize = 24)
    ax.set_xlabel('RA', fontsize = 20)
    ax.set_ylabel('Dec', fontsize = 20)
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    ax.tick_params(axis='both', labelsize=18)
    
    # Plot and save

    return ax






        
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
    axes[0,0].set_xlabel(r'$V$ Johnson-Cousin')
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
    axes[0,1].set_xlabel(r'$V$ Johnson-Cousin')
    axes[0,1].set_ylabel('Number of stars')
    axes[0,1].tick_params(axis='x', which='minor', bottom=True, top=False)
    axes[0,1].tick_params(axis='x', which='major', bottom=True, top=False)
    axes[0,1].tick_params(axis='y', which='minor', left=False, right=False)
    axes[0,1].tick_params(axis='y', which='major', left=True, right=False)
    axes[0,1].grid(axis='y', color='gray', alpha=0.3)

    # Prepare bins and plot number distribution of contaminants per target

    numbinCon  = 1 + int(np.max(numConPerTar)/100)
    binsizeNum = int((np.max(numConPerTar) - 0) / numbinCon) + 2
    binlistNum = np.linspace(-0.5, np.max(numConPerTar)+0.5, binsizeNum)  # -0.5 because num x-axis
    axes[1,0].hist(numConPerTar, binlistNum, facecolor='g', edgecolor='g', fill=True, log=True, alpha=0.3)
    axes[1,0].yaxis.set_major_formatter(ScalarFormatter())
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

    # Finito!

    return axes









def plotYawPitchRollTimeSeries(fig, time, signals, units, title=False, ylim=False):
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

    try:
        numData = len(signals)
    except:
        numData = 1

    # Handle yaxis limits

    if ylim is False:
        ylim = 0.5*np.max(np.abs(signals))

    # Adjust linewidth after data

    if len(time) < 1e3: lw = 1
    else: lw = 0.5

    # Make plot

    labels = ['Yaw', 'Pitch', 'Roll']
    colors = ['royalblue', 'lightseagreen', 'limegreen']

    for plot in range(numData):

        axes = fig.add_subplot(numData, 1, plot+1)

        # Make sure that time series is residuals around zero

        signals[plot] -= np.median(signals[plot])

        # Plot timeseries

        axes.plot(time, signals[plot], '-', c=colors[plot], lw=lw)

        # Add root-mean-square lines

        rms = np.sqrt(np.mean(signals[plot]**2))
        axes.axhline(+rms, c='k', ls='--', lw=0.7, label='RMS = {0:.3f} {1}'.format(rms, units[1]))
        axes.axhline(-rms, c='k', ls='--', lw=0.7)
        axes.legend(loc='upper right')

        # Latter settings

        axes.set_ylabel('{0} [{1}]'.format(labels[plot], units[1]))
        axes.set_xlim(np.min(time), np.max(time))
        axes.set_ylim(-ylim, +ylim)

        # Remove tick labels on x axis except for last plot

        if plot < numData-1:
            axes.tick_params(labelbottom=False)
            #axes.set_ylim(-lim, +lim)
        else:
            axes.set_xlabel('Time [{0}]'.format(units[0]))

        # Title

        if plot == 0: axes.set_title(title, fontsize=fs)

    # Adjust layout

    plt.tight_layout()
    plt.subplots_adjust(hspace = .001)

    # Finito!

    return axes









def plotYawPitchRollPSD(fig, time, signals, scale=1e-6, carbox=144, title=False, labels=False, xmin=False, ylim=False, misreq=False):
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
    carbox : int (optional)
        Length of median carbox filter. Default is 3600s/25s = 144
    title : str (optional)
        Title for plot
    labels : list-str (optinal)
        List of string labels where the first is the xlabel and the rest is ylabels
    xmin : float (optional)
        Limit for x min. The x max limit is the Nyquist frequency
    ylim : list-float (optional)
        List of y min and max limit ["y-min", "y-max"]
    misreq : bool (optional)
        If "True" the mission requirements for the AOSC will be plotted alonside the data.

    Return
    ------
    Plot or/and saved plot to PNG.
    """

    # Number of data sets

    numData = len(signals)

    # Find time step

    sampling = (time[1]-time[0]) * scale

    # Make plot

    labels = ['Yaw', 'Pitch', 'Roll']
    #colors = ['tomato', 'darkorange', 'gold']
    colors = ['royalblue', 'lightseagreen', 'limegreen']

    for plot in range(numData):

        # Create axes objects

        axes = fig.add_subplot(numData, 1, plot+1)

        # Find PSD and median filter

        freq, PSD = ns.powerDensityFFT(signals[plot], sampling)
        PSD_med   = median_filter(PSD, carbox)
        perhour   = int(carbox*sampling/3600.)

        # Plot results

        axes.plot(freq, PSD,     '-', c=colors[plot], lw=lw, label=labels[plot])
        axes.plot(freq, PSD_med, 'k-', lw=lw+1, label='Median filter')#label='{0}h median '.format(perhour))

        # Plot mission requirements (from the red book)

        if misreq:
            axes.plot([3e-6*scale, 20e-6*scale], [21.4*scale, 0.23*scale], c='k', linestyle='--', lw=1, label='MPE requirement')
            axes.plot([20e-6*scale, 4e-2*scale], [0.23*scale, 0.23*scale], c='k', linestyle='--', lw=1)

        # Log scaling

        axes.set_xscale("log")
        axes.set_yscale("log")

        # Latter settings

        if plot == 1:
            axes.set_ylabel(r'Amplitude [arcsec$^2$ Hz$^{-1}$]')

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

        axes.legend(loc='best')

        # Set title

        if title is not False and plot == 0: axes.set_title(title, fontsize=fs)

    # Remaining

    plt.xlabel(r'Frequency [Hz]')
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

    lim    = np.max(np.abs(signals))
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
        lim = np.max(np.abs(signals))
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









def plotPSD(fig, freq, psd, carbox=144, units=False, labels=False, colors=False, title=False,
            xlim=False, ylim=False, linewidth=False, misreq=False):
    """
    This function plots the Power Spectral Density (PSD). Alongside the data a median
    filter is plotted with a default carbox length of 144 time points, corresponding to
    1 hour precision if the time series has a cadence of 25 seconds.

    Parameters
    ----------
    freq : narray
        Frequency points [Hz, mHz, or mizroHz]
    psd : narray, list-narray
        Either single signal array or a list of signal arrays
    carbox : int (optional)
        Length of median carbox filter. Default is 3600s/25s = 144. Also False to ignore.
    title : str (optional)
        Title for plot
    labels : list-str (optinal)
        List of string labels where the first is the xlabel and the rest is ylabels
    xmin : float (optional)
        Limit for x min. The x max limit is the Nyquist frequency
    ylim : list-float (optional)
        List of y min and max limit ["y-min", "y-max"]

    Return
    ------
    Plot or/and saved plot to PNG.
    """

    # Handle the number of input data sets

    if type(psd) == list:
        numData = len(psd)
    else:
        numData = 1
        freq = [freq]
        psd = [psd]

    # Handle axes units

    if units is False:
        units = ['$\mu$Hz', 'ppm']
        scale = 1e6
    else:
        scale = 1

    # Handle colors

    if colors is False:
        colors = ['tomato', 'darkorange', 'gold']
        if numData > 3:
            colors = cm.rainbow(np.linspace(0, 1, numData))

    # Handle linewidths

    if linewidth is False:
        lw = 1
    else:
        lw = linewidth

    # Create axes objects

    axes = fig.add_subplot()

    # Allow plotting multiple PSDs in consecutive subplots

    for plot in np.arange(numData):

        # Plot results

        if labels is False:
            plt.plot(freq[plot], psd[plot], '-', c=colors[plot], lw=lw)
        else:
            plt.plot(freq[plot], psd[plot], '-', c=colors[plot], lw=lw, label=labels[plot])

        # Plot median filter if requested

        if carbox:
            perhour = carbox*25/3600
            PSD_med = median_filter(psd[plot], carbox)
            plt.plot(freq[plot], PSD_med, 'k-', lw=lw+1, label='{0}h median '.format(perhour))

    # Plot mission requirements (from the red book)

    if misreq:
        plt.plot([3e-6*scale, 20e-6*scale], [21.4*scale, 0.23*scale], c='k', linestyle='--', lw=1, label='MPE requirement')
        plt.plot([20e-6*scale, 4e-2*scale], [0.23*scale, 0.23*scale], c='k', linestyle='--', lw=1)

    # Log scaling

    plt.xscale("log")
    plt.yscale("log")

    # Latter settings

    plt.ylabel(r'PSD [{}$^2$ {}'.format(units[1], units[0])+'$^{-1}$]')

    # Set x-min limit

    if xlim is not False:
        plt.xlim(xlim[0], xlim[1])

    # Set y limits

    if ylim is not False:
        plt.ylim(ylim[0], ylim[1])

    # Set title

    if title is not False and plot == 0:
        plt.title(title, fontsize=fs)

    if labels is not False or misreq is True:
        plt.legend(loc='best')

    # Remaining

    plt.xlabel(r'Frequency [{}]'.format(units[0]))
    plt.tight_layout()
    plt.subplots_adjust(hspace = .001)
    plt.grid()

    # Finito!

    return axes








def plotPhotometry(fig, outputFile, medfilt=144, fluxInput=False, NSR=False, COB=False, title=False):
    """
    PURPOSE: Plot the photometric time series of single target from a HDF5 file.

    PARAMETERS
    ----------
    outputFile : str
        File name of HDF5 file containing photometry
    medfilt : int
        Number of time points for overlaid median filter.
        Default is 1h cadence: 3600s/25s=144 time points.
    NSR : bool
        If the Noise-to-Signal Ratio (NSR) should be plotted alonside the time series.
        Default is False.
    COB : bool
        If the Center-Of-Brightness (COB) should be plotted alonside the time series.
        Default is False.
    title : str
        String of title plot. Default is False.

    RETURN
    ------
    axes : object
        Axes matplotlib.pyplot handle object to be modified by the user
    """

    # Load photometry class

    f = SimFile(outputFile)

    # Fetch photometry

    lc = f.getPhotometry(1)
    time     = lc[0] / c.day        # [days]
    flux_in  = ut.normalize(lc[1])  # [ppm]
    flux_out = ut.normalize(lc[2])  # [ppm]

    # Compute median carbox filter of output signal

    flux_med = median_filter(flux_out, medfilt)  # [ppm]

    # Fetch mask updates

    mask = f.getPhotometricMask(1)
    maskupdates = (mask[2] * 25.) / c.day   # [days]
    maskNSR     = mask[4] * 1e6             # [ppm]

    # Create figure and subplots

    plt.subplots_adjust(wspace=0.15, hspace=0.20)

    # Plot photometric time series

    ax0 = fig.add_subplot(2,1,1)
    ax0.plot(time, flux_out, 'k.', markersize=1, alpha=0.2, label='Raw flux')
    ax0.plot(time, flux_med, 'g-', label='Median per hour')
    if fluxInput is not False: ax0.plot(time, flux_in, '-', c='royalblue', markersize=2, label='Model')
    axes_maskupdates(ax0, time, maskupdates)
    ax0.legend(loc='lower right', fancybox=True, ncol=2)
    ax0.set_ylabel('Norm. Flux [ppm]')

    # Plot NSR

    if NSR is not False:

        ax1 = fig.add_subplot(2,1,2)
        ax1.plot(maskupdates, maskNSR, 'k--', alpha=0.5)
        ax1.plot(maskupdates, maskNSR, 'm*', label='Mask update NSR')
        ax1.legend(loc='upper right',  fancybox=True)
        ax1.set_xlabel('Time [days]')
        ax1.set_ylabel(r'NSR [ppm h$^{-1}$]')
        ax1.set_xlim(axes_minmax(x=time))

    if COB is not False:

        # Fetch pixel coordinates

        rowPix, colPix = f.getStarPositions(1)

        # Plot row pixel on left y axis

        ax1 = fig.add_subplot(2,1,2)
        ax1.plot(time, rowPix, '-', c='darkcyan', label='Row pixel')
        ax1.set_ylabel('Row coordinate [pixel]')
        ax1.set_xlabel('Time [days]')
        ax1.set_title('Star position')
        ax1.set_xlim(axes_minmax(x=time))
        ax1.legend(loc='upper left')
        
        # Plot column pixel on right y axis

        ax2 = ax1.twinx()
        ax2.plot(time, colPix, '-', c='hotpink', label='Col pixel')
        ax2.set_ylabel('Column coordinate [pixel]')
        ax2.set_xlim(axes_minmax(x=time))
        ax2.legend(loc='lower right')
        
    # Labels

    if title is not False:
        fig.text(0.5, 0.9, title, ha='center', fontsize=fs)
    if NSR is False and COB is False:
        ax0.set_xlabel('Time [days]')
        
    # Limits

    flux_min = np.min([flux_in, flux_out])
    ax0.set_ylim(flux_min, np.abs(flux_min))
    ax0.set_xlim(axes_minmax(x=time))

    # Finito!

    plt.show()












def plotMultiCameraAndQuarterPhotometry(fig, outputFiles, medfilt=144, title=None):
    """
    PURPOSE: 

    PARAMETERS
    ----------
    filename : str
        File name of HDF5 file containing photometry

    RETURN
    ------
    Plot or/and saved plot to PNG.

    FIXME this function do not work currently. Perhaps this should be moved to photometry file
    as a bigger library of function to plot multi-quarter and multi-camera time series
    """

    # Fetch information about the observation

    numFiles = len(outputFiles)
    aa = 0.2  # Alpha channel of data

    ax = fig.add_subplot(1,1,1)

    cameras = {}
    maskupdates = []

    for i in range(numFiles):

        # Load photometry class

        f = SimFile(outputFiles[i])

        quarter_i = int(outputFiles[i][-6])
        camera_i  = outputFiles[i][-11:-8]

        #print(quarter_i)
        # Fetch and combine time series for each quarter

        # Fetch photometry

        lc = f.getPhotometry(starID=1, quarterNo=quarter_i)
        time = lc[0] / c.day        # [days]
        flux = ut.normalize(lc[2])  # [ppm]
        flux_med = median_filter(flux, medfilt)  # [ppm]

        ax.plot(time, flux, 'k.', markersize=1, alpha=0.2, label='Raw flux')
        ax.plot(time, flux_med, 'g-', label='Median per hour')

        # Fetch mask updates

        mask = f.getPhotometricMask(1)
        maskupdates.append((mask[2] * 25.) / c.day)  # [days]

            # Save time series to dict

        #cameras['Ncam{}'.format(camera_i)] = [time, flux]


    # Compute median carbox filter of output signal

    maskupdates = np.arange(0, 7*8*maskupdates[0][1], maskupdates[0][1])

    axes_maskupdates(ax, time, maskupdates)
    ax.legend(loc='lower right', fancybox=True, ncol=2)
    ax.set_ylabel('Norm. Flux [ppm]')
    plt.xlabel('Time [days]', fontsize=fs-3)
    plt.ylabel('Relative flux [ppm]', fontsize=fs-3)

    # Settings

    #plt.xlim(axes_minmax(x=time))
    #plt.ylim(-n*df, df)

    # That's it!

    plt.show()











def plotPhotometryComparison(fig, filenames, medfilt=None, title=None):
    """
    PURPOSE: 

    PARAMETERS
    ----------
    filename : str
        File name of HDF5 file containing photometry

    RETURN
    ------
    Plot or/and saved plot to PNG.

    FIXME this function do not work currently. Perhaps this should be moved to photometry file
    as a bigger library of function to plot multi-quarter and multi-camera time series
    """

    # User defined labels

    #title  = 'Drift test: 10.0 mag; All CCD effects ON; Jitter OFF'
    #labels = ['0.0', '0.5', '1.0', '1.5', '2.0', '2.5']

    #title  = 'Drift test: 10.0 mag; All CCD effects ON; Jitter of 0.04 RMS'
    #labels = ['0.1', '2.0', '3.0']

    # title  = 'Jitter test: All CCD effects on & drift off'
    # labels = ['0.00', '0.01', '0.02', '0.03', '0.04', '0.05', '0.06']

    #title  = 'Test of flux decrease: 10.0 mag; WASP-33b hot-Jupiter'
    #labels = ['TE con. + Drift/jitter OFF', 'TE reg. + Drift/jitter OFF', 'TE con. + 2.0/0.04 RMS', 'TE req. + 2.0/0.04 RMS']

    #title  = 'Magnitude test: All CCD effects ON; Drift 2.0 RMS; Jitter 0.04 RMS'
    #labels = ['8.0 mag', '9 mag', '10 mag']

    #title  = 'Test of exoplanet input model: 10.0 mag; Sun + hot-Jupiter'
    #labels = ['Disabled: Drift, Jitter, CTI, and TE']

    #title  = 'Test of exoplanet input model: 10.0 mag; Constant flux'
    #labels = ['Disabled: Drift, Jitter, CTI, and TE', 'Disabled: Drift, Jitter, CTI', 'Disabled: Drift, Jitter']

    #title  = 'Test Photometry: Constant flux; Jitter 0.04 arcsec RMS; TED yaw and pitch = 15 arcsec/quarter, roll from Prime'
    #labels = ['10.0 mag', '11.0 mag']

    # Fetch information about the observation

    n = len(filenames)

    method = 'median'
    maskupdate = 14
    df = 4e4
    aa = 0.2

    # Colvolution filter

    if medfilt is None: medfilt = 144

    # Plot each time series with an offset

    for i in range(n):

        photometryClass = PhotometricFile(filenames[i])
        signals = photometryClass.getPhotometricTimeSeries(1)

        if i == 0: time = signals[0]/c.day  # [days]

        flux = normalize(signals[2]) - df*(i)
        plt.plot(time, flux, 'o', c='k', markersize=1, alpha=aa)

        if method == 'model':

            flux_input = normalize(signals[1]) - df*(i)
            plt.plot(time, flux_input, '-', c=cb[i+1], label=labels[i])
            plt.axhline(y=-i*df, c='gray', linestyle='--', linewidth=1)

        if method == 'median':

            flux_med = median_filter(flux, medfilt)
            plt.plot(time, flux_med, '-', c=cb[i+1], markersize=1, label=labels[i])
            plt.axhline(y=np.median(flux_med[:1000]), c='gray',linestyle='--', linewidth=1)

    # Plot Quarters

    quarters = np.arange(0, time[-1]/(24*3600.), 120)
    for Q in quarters:
        if Q == 0:
            plt.axvline(x=Q, c='darkgray', linestyle='-.', linewidth=1, label='Quarter marks')
        else:
            plt.axvline(x=Q, c='darkgray', linestyle='-.', linewidth=1)

    # Plot mask update occurance

    updates = np.arange(0, 90, maskupdate)
    for update in updates:
        if update == 0:
            plt.axvline(x=update, c='k', linestyle=':', linewidth=1, label='Mask updates')
        else:
            plt.axvline(x=update, c='k', linestyle=':', linewidth=1)

    # Labels

    if title is not None: fig.text(0.5, 0.9, title, ha='center', fontsize=fs)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 0.995), fancybox=True, ncol=n+2, fontsize=fs-3)
    plt.xlabel('Time [days]', fontsize=fs-3)
    plt.ylabel('Relative flux [ppm]', fontsize=fs-3)

    # Settings

    plt.xlim(axes_minmax(x=time))
    plt.ylim(-n*df, df)

    # That's it!

    plt.show()










def plotSubfieldAnimation(fig, filename, numImages=False, outputFileName=False, clipPercentile=8.0,
                          useTitle=True, colorMap="hot", showGrid=False, showStarPositions=False,
                          showPointLikeGhostPositions=False, minVmag=None, maxVmag=None,
                          showStarIDs=False, showMaskOfStarID=None, skipNimages=None):
    """
    Create and plot an animation of a set of imagettes.

    Parameters
    ----------
    filename : str

    Return
    ------
    axes : object
        Axes matplotlib.pyplot handle object to be modified by the user

    EXAMPLES
    --------
    fig = plt.figure(fig=(10,10))
    TODO add photometric mask to animation
    """

    # Fetch file with simulated subfields

    f = h5py.File(filename, "r")

    # Fetch the image names
    imgNames   = list(f['Images'].keys())
    #imgNumbers = list(range(1, len(imgNames)+1))
    imgNumbers = list(range(len(imgNames)))
    N          = len(imgNames)

    # Skip N images to make animation faster
    if skipNimages is not None:
        imgNames   = imgNames[0::skipNimages]
        imgNumbers = imgNumbers[0::skipNimages]
        
    # Plot the image. Note that pixel coordinates start at the left bottom side of each pixel.

    ims = []
    fig, axis = plt.subplots(1,1)
    for imgNumber, imgName in zip(imgNumbers, imgNames):

        # Fetch each pixel image
        image = f["Images/{0}".format(imgName)]

        # Fetch image dimention first time only
        if not ims:
            Nrows, Ncols = image.shape

        # Make image plot
        imagePlot = axis.imshow(image, cmap=colorMap, interpolation="nearest",
                                origin='lower', extent=[0,Nrows,0,Ncols], animated=True)

        # The large dynamic range of the pixel values often results in images where only
        # the brightest stars are visible. To improve the contrast, clip the color mapping.
        imagePlot.set_clim(np.percentile(image, clipPercentile), np.percentile(image, 100-clipPercentile))

        # OVERPLOT STAR POSITIONS

        if showStarPositions:
            # Extract the arays from HDF5 file
            exposureGroupName = "Exposure{0:06d}".format(imgNumber)
            dataset = f["StarPositions"][exposureGroupName]["starID"]
            ID = np.zeros(dataset.shape, dataset.dtype)
            dataset.read_direct(ID)
            dataset = f["StarPositions"][exposureGroupName]["rowPix"]
            row = np.zeros(dataset.shape, dataset.dtype)
            dataset.read_direct(row)
            dataset = f["StarPositions"][exposureGroupName]["colPix"]
            col = np.zeros(dataset.shape, dataset.dtype)
            dataset.read_direct(col)
            dataset = f["StarPositions"][exposureGroupName]["flux"]
            flux = np.zeros(dataset.shape, dataset.dtype)
            dataset.read_direct(flux)

            # Allow differentiating between a (PIC) target and its contaminants
            if showStarPositions == 'PIC':
                tarMarkerSize = 200
                mag = -2.5*np.log10(flux)
                coor_tar = axis.scatter(col[0], row[0], s=tarMarkerSize, marker='o', c='lime',
                                        edgecolor='k', linewidth=1, zorder=4)
                if len(col) > 1:
                    conMarkerSize = (tarMarkerSize /
                                     (mag[1:] - mag[0]*np.ones(len(col)-1))).astype(int)
                    coor_con = axis.scatter(col[1:], row[1:], s=conMarkerSize, marker='o', c='gold',
                                            edgecolor='k', linewidth=1, zorder=4)
            # Or hightligth all stars the same
            else:
                axis.scatter(col, row, marker='x', c='g')
            if showStarIDs:
                for k in range(len(ID)):
                    label = "{0}".format(ID[k])
                    axis.annotate(label, (col[k], row[k]), fontsize='small', fontweight='extra bold', color="black")
                    
        # Ensure that the axis limits are properly set
        axis.set_xlim(0, Ncols)
        axis.set_ylim(0, Nrows)

        # User defined title-string
        if isinstance(useTitle, str):
            plt.title(useTitle)

        # By default, matplotlib only shows the (x,y) coordinates of each pixel but not
        # the pixel value itself. Change this by redefining the axis.format_coord
        def format_coord(x, y):
            col = int(x)
            row = int(y)
            if col >= 0 and col < Ncols and row >= 0 and row < Nrows:
                z = image[row,col]
                return "x={:.1f}, y={:.1f}, z={:.1f}".format(x, y, z)
            else:
                return "x={:.1f}, y={:.1f}".format(x, y)
        axis.format_coord = format_coord

        # Show all ticks for smaller subfields or otherwise 10
        if Ncols < 10 and Nrows < 10:
            plt.xticks(np.arange(0, Nrows+1))
            plt.yticks(np.arange(0, Ncols+1))
        else:
            plt.xticks(np.arange(0, Nrows, 10))
            plt.yticks(np.arange(0, Ncols, 10))
            plt.tight_layout()
        # Overplot rectangles over those pixels that are part of the mask
        # Note: imshow reverses rows and columns

        # if showMaskOfStarID is not None:
        #     rowIndices, colIndices, exposureNr = getPhotometricMask(showMaskOfStarID, imageNr)
        #     for k in range(len(rowIndices)):
        #         rect = patches.Rectangle((colIndices[k], rowIndices[k]), 1, 1, linewidth=2.0, edgecolor='b', facecolor='none')
        #         axis.add_patch(rect)

        # If requiered, overplot a gray semi-transparent grid
        # Note: this is only meaningsful for smaller imagettes

        if showGrid is True:
            axis.grid(c='gray', ls='-', alpha=0.3)

        # Add x and y axis labels
        plt.xlabel('x [pixel]')
        plt.ylabel('y [pixel]')
            
        # Append images to list
        ims.append([imagePlot, coor_tar, coor_con])

        # Compile to bash
        ut.compilation(imgNumber, N, 'Done')
    print; print('')

    # CREATE ANIMATION
    ani = animation.ArtistAnimation(fig, ims, interval=100, blit=True, repeat_delay=0)
    
    # Save animation (fps=50 and dpi=100 seems like good settings)
    print('Creating GIF animation..')
    if outputFileName is not False:
        ani.save(f'{outputFileName}.gif', fps=50, dpi=100)

    # Show animation
    plt.draw()
    plt.plot()

    # Finito!
    return

