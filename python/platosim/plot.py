#!/usr/bin/env python3

import os
import h5py
import numpy as np

from matplotlib.pyplot import cm
from matplotlib import patches
from matplotlib.path import Path
from matplotlib.ticker import MaxNLocator, ScalarFormatter, LogLocator
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from scipy import constants as c
from scipy.ndimage import median_filter

import astropy.units as u
from astropy.coordinates import SkyCoord

import ligo.skymap.plot
import shapely.geometry as sg
import descartes

# PlatoSim modules

import platosim.referenceFrames as rf
import platosim.utilities       as ut
import platosim.noise           as ns
from platosim.matplotlibrc import setup
setup()

# Constants

aa = 0.5  # Alpha
ms = 3    # MarkerSize
lw = 2    # LineWidth


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



def discrete_colorbar(cbins, cmap="coolwarm"):
    """
    """
    import matplotlib as mpl
    from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
    
    # define the colormap
    cmap = plt.get_cmap(cmap)
    
    # Extract all colors from the the choosen colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]

    # Create the new map
    cmap = LinearSegmentedColormap.from_list('custom', cmaplist, cmap.N)

    # define the bins and normalize
    norm = mpl.colors.BoundaryNorm(cbins, cmap.N)

    return norm




def slider(imagePlot, images, Nimg, label="Image number"):

    # Function to update slider
    def update_image(n=0):
        image = images[n]
        imagePlot.set_data(image)
        fig.canvas.draw()


    # Show the image
    slider = IntSlider(0, 0, 10, 1, layout=Layout(width='500px'))

    from ipywidgets import interact
    interact(update_image, n=(0,Nimg-1), x=slider)

    # Make a horizontal slider to control the frequency.
    #axfreq = fig.add_axes([0.25, 0.1, 0.65, 0.03])
    #slider = Slider(ax=axfreq, label=label, valmin=vmin, valmax=vmax, valinit=val)


#==============================================================#
#                         PROJECTIONS                          #
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











def drawStarInCCDfocalPlane(fig, sim, xCCD, yCCD, refCcdCode, refGroup, raPlatform, decPlatform, tiltAngle, azimuthAngle, solarPanelOrientation):
    """
    PURPOSE:  Draw a star given by the CCD pixel coordinates in the CCD focal plane.

              This function plot the CCD pixel position in the CCD focal plane, together
              with the CCD- and Camera footprints. Given the reference CCD- and camera group,
              the star location is indicated by black arrows sprung from the origo of the 
              parrent CCD. To easy the eye (of this rather complicated plot) the star is
              highligthed by the parrent CCD color and the dashed line is connecting the
              star to the optical axis of the camera-group. 

    INPUT:    fig:                    matplotlib figure object (e.g.: fig = fig.plt.figure(figsize=(10,10)))
              sim:                    instance of simulation class (see simulation.py)
              xCCD:                   column  pixel coordinate [pixel]
              yCCD:                   row pixel coordinate [pixel]
              refCcdCode:             Reference of CCD (1, 2, 3, 4)
              refGroup:               Reference of camera group (1, 2, 3, 4)
              raPlatform:             Right Ascension of platform pointing [deg]
              decPlatform:            Declination of platform pointing [deg]
              titlAngle:              Tilt angle of camera [deg]
              azimuthAngle:           Azimuth angle of camera [deg] 
              solarPanelOrientation:  Orientation of solar panel [deg]
              
    OUTPUT:   Plot only, no axes object.

    NOTE:     This plot neglects pointing- and camera alignment errors.
    """

    # Parameter used for the plot
    
    numGroups     = 4
    numCorners    = 4
    offset        = 4
    colors        = ['b', 'r', 'g', 'orange']
    ccdCodes      = ["1", "2", "3", "4"]
    tiltAngles    = sim['CameraGroups/TiltAngle'][:numGroups]           # [deg]
    azimuthAngles = sim['CameraGroups/AzimuthAngle'][:numGroups]        # [deg]
    fovDegrees    = sim['CCD/RelativeTransmissivity/RadiusFOV']         # [deg]
    focalLength   = sim['Camera/FocalLength/ConstantValue'] * 1e3       # [mm]
    pixelSize     = sim['CCD/PixelSize']                                # [micron]
    plateScale    = sim['Camera/PlateScale'] * pixelSize                # [arcsec]

    # Find actual FOV in pixel and mm
    
    fovPixels  = fovDegrees / plateScale * c.degree / c.arcsec
    fovMm      = focalLength * np.tan(np.radians(fovDegrees))

    def mm2pixels(distanceMm):
        """
        Conversion from millimeters to pixels.
        :param distanceMm: Distance [mm].
        :return distancePixels: Distance [pixels].
        """
        distancePixels = np.degrees(np.arctan(distanceMm / focalLength)) / plateScale * c.degree / c.arcsec
        return distancePixels

    sign = lambda x: (1, -1)[x < 0]

    xFP = np.array([])
    yFP = np.array([])
    # Position of the Sun
    raSun, decSun = rf.sunSkyCoordinatesAwayfromPlatformPointing(np.radians(raPlatform),
                                                                 np.radians(decPlatform),
                                                                 np.radians(solarPanelOrientation))

    # Telescope pointing w.r.t. platform pointing
    
    raTelescope = []
    decTelescope = []
    for group in range(0, numGroups):
        
        # Telescope pointing (absolute) [radians]
        
        ra, dec = rf.platformToTelescopePointingCoordinates(np.radians(raPlatform), np.radians(decPlatform),
                                                            raSun, decSun,
                                                            np.radians(azimuthAngles[group]),
                                                            np.radians(tiltAngles[group]))

        # Telescope pointing w.r.t. platform pointing
        
        raTelescope.append(np.degrees(ra) - raPlatform)     # [degrees]
        decTelescope.append(np.degrees(dec) - decPlatform)  # [degrees]

    meanDist = (np.mean(np.absolute(raTelescope)) + np.mean(np.absolute(decTelescope))) / 2.0
    for group in range(numGroups):
        raTelescope[group]  = sign(raTelescope[group]) * meanDist
        decTelescope[group] = sign(decTelescope[group]) * meanDist

    # Calculate the coordinates of the telescope pointings (for the 4 telescope groups) 
    # in the focal-plane reference frame of group 1
    
    for group in range(numGroups):
        xFP1, yFP1 = rf.skyToFocalPlaneCoordinates(np.radians(raTelescope[group] + raPlatform), 
                                                   np.radians(decTelescope[group] + decPlatform), 
                                                   np.radians(raPlatform), np.radians(decPlatform),
                                                   np.radians(solarPanelOrientation),
                                                   np.radians(tiltAngles[0]), 
                                                   np.radians(azimuthAngles[0]), 
                                                   0, focalLength)
        xFP = np.append(xFP, xFP1)
        yFP = np.append(yFP, yFP1)

    # Make sure the average of the telescope pointings of the 4 telescope groups is at the origin of the reference frame
    
    xAvg = np.mean(xFP)
    yAvg = np.mean(yFP)
    xFP -= xAvg
    yFP -= yAvg

    meanDist = (np.mean(np.absolute(xFP)) + np.mean(np.absolute(yFP))) / 2.0

    for group in range(numGroups):
        xFP[group] = sign(xFP[group]) * meanDist
        yFP[group] = sign(yFP[group]) * meanDist

    xPixels = np.copy(xFP)
    yPixels = np.copy(yFP)

    for group in range(numGroups):
        xPixels[group] = mm2pixels(xPixels[group])     # [mm] -> [pixels]
        yPixels[group] = mm2pixels(yPixels[group])     # [mm] -> [pixels]

    index = 0

    cornersX, cornersY = rf.computeCCDcornersInFocalPlane(refCcdCode, pixelSize)
    offsetX = mm2pixels(cornersX[index]) + xPixels[refGroup - 1]
    offsetY = mm2pixels(cornersY[index]) + yPixels[refGroup - 1]

    # Correct input pixel coordinates to match orientation of CCD origin
    if refCcdCode == '1': xCCD, yCCD = +xCCD, +yCCD
    if refCcdCode == '2': xCCD, yCCD = -yCCD, +xCCD 
    if refCcdCode == '3': xCCD, yCCD = -xCCD, -yCCD
    if refCcdCode == '4': xCCD, yCCD = +yCCD, -xCCD
    
    # START PLOT
    
    ax = fig.add_subplot(111)

    # Add grid

    ax.grid(True)

    # Gray shade N-CAM visibility

    circles = []
    for group in range(numGroups):
        circles.append(sg.Point(-(xPixels[group] - offsetX), -(yPixels[group] - offsetY)).buffer(fovPixels))
    for index in range(numCorners):
        one = circles[index].intersection(circles[index])
        ax.add_patch(descartes.PolygonPatch(one, fc='gray', ec='none', alpha=0.2))
        two = circles[index].intersection(circles[(index + 1) % numCorners])
        ax.add_patch(descartes.PolygonPatch(two, fc='gray', ec='none', alpha=0.1))
        three = circles[index].intersection(circles[(index + 1) % numCorners]).intersection(circles[(index + 2) % numCorners])
        ax.add_patch(descartes.PolygonPatch(three, fc='gray', ec='none', alpha=0.05))
    four = circles[0].intersection(circles[1]).intersection(circles[2]).intersection(circles[3])
    ax.add_patch(descartes.PolygonPatch(four, fc='gray', ec='none', alpha=0.03))

    # Plot CCD footprint ontop

    for ccdCode in ccdCodes:
        cornersX, cornersY = rf.computeCCDcornersInFocalPlane(ccdCode, pixelSize)
        for corner in range(numCorners):
            cornersX[corner] = mm2pixels(cornersX[corner])
            cornersY[corner] = mm2pixels(cornersY[corner])
        cornersX = np.append(cornersX, cornersX[0])     # [mm]
        cornersY = np.append(cornersY, cornersY[0])     # [mm]
        for group in range(numGroups):
            arrayX = cornersX + xPixels[group] - offsetX
            arrayY = cornersY + yPixels[group] - offsetY
            for index in range(5):
                arrayX[index] = (arrayX[index])
                arrayY[index] = (arrayY[index])
            if group != refGroup-1:
                plt.plot(-arrayX, -arrayY, color=colors[group], zorder=1)
            else:
                plt.plot(-arrayX, -arrayY, color=colors[group], zorder=2)
            off = 4000
            if group == refGroup-1 and int(ccdCode) == int(refCcdCode):
                ax.text(-(xPixels[group] - offsetX + sign(cornersX[3]) * off),
                        -(yPixels[group] - offsetY + sign(cornersY[3]) * off),
                        ccdCode, fontsize=20, color=colors[group], weight='bold')
            else:
                ax.text(-(xPixels[group] - offsetX + sign(cornersX[3]) * off),
                        -(yPixels[group] - offsetY + sign(cornersY[3]) * off),
                        ccdCode, fontsize=15, color=colors[group])

    # Plot center of group

    for group in range(numGroups):
        plt.plot([-(xPixels[group] - offsetX)], [-(yPixels[group] - offsetY)], mfc=colors[group], marker="o", mec='k', ms=10)
        circ = plt.Circle(((-(xPixels[group] - offsetX)), (-(yPixels[group] - offsetY))),
                          radius=fovPixels, color="none", linewidth=1, label="Group " + str(group + 1))
        ax.add_patch(circ)
        circ.set_edgecolor(colors[group])
        circ.set_facecolor("none")

    # Plot arrow

    xHead, yHead = -200, -200
    if xCCD < 0: xHead *= -1
    elif xCCD < 200: xHead = abs(xCCD)
    if yCCD < 0: yHead *= -1
    elif yCCD < 200: yHead = abs(yCCD)
    ax.arrow(0, 0, xCCD+xHead, 0, head_width=150, head_length=abs(xHead), fc='k', ec='k', linewidth=3, zorder=5)
    ax.arrow(0, 0, 0, yCCD+yHead, head_width=150, head_length=abs(xHead), fc='k', ec='k', linewidth=3, zorder=5)
    plt.plot([0], [0], "ko", ms=10, label="CCD origin", zorder=6)
        
    # Plot target star

    plt.plot([-(xPixels[refGroup-1] - offsetX), xCCD], [-(yPixels[refGroup-1] - offsetY), yCCD], 'k--', zorder=6)
    plt.plot(xCCD, yCCD, '*', mfc=colors[refGroup-1], mec='k', ms=19, label='Target star', zorder=7)
    
    # Settings

    ax.set_title("CCD " + refCcdCode + " in Group " + str(refGroup), fontsize=20)
    plt.legend(prop={"size": 12}, bbox_to_anchor=(1.0, 1.0))
    ax.set_xlabel("Column [pixel]", fontsize=15)
    ax.set_ylabel("Row [pixel]", fontsize=15)
    ax.set_aspect('equal', 'box')
    plt.tight_layout()

    # Finito!
    
    plt.show()













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
    fs = 20
    plt.title('Aitoff projection in Galactic coordinates', fontsize=fs+2, y=1.02)
    fig, ax = fig
    if len(raStars) <= 1e2: ms = 3.
    if len(raStars) >= 1e2 and len(raStars) < 1e3: ms = 1.3
    if len(raStars) >= 1e3 and len(raStars) < 1e5: ms = 1.
    if len(raStars) >= 1e5: ms = 0.1

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
        cbarax = fig.add_axes([0.805, 0.2, 0.02, 0.57])
        cbar = plt.colorbar(im, orientation='vertical', cax=cbarax, extend='both')
        cbar.set_label(r'PLATO passband, $P$', fontsize=fs)
        cbar.ax.tick_params(labelsize=fs)
    else:
        cbarax = fig.add_axes([0.25, 0.08, 0.525, 0.03])
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

    # Add axes labels
    axes.set_xlabel(r'Longitude, $l$ [deg]', fontsize=fs)
    axes.set_ylabel(r'Latitude, $b$ [deg]', fontsize=fs)

    # Set grid and remore outer ticks (if set by default)
    axes.grid(True, alpha=0.3)
    ax.axis('off')
    plt.draw()
    plt.tight_layout()

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







def plotPlatoFOV(pointingField, raStars, decStars, magStars=None,
                 nCamVis=None, skymap=None, title=None, fs=20):
    """
    Funtion to plot 

    Parameters
    ----------

    Return
    ------
    fig : object
        Axes matplotlib.pyplot handle object to be modified by the user
        Use fig..savefig('<plot.png>', bbox_inches='tight', dpi=200)
    """

    # Select field
    indir = os.getenv('PLATO_PROJECT_HOME') + '/python/platosim/picsim/'
    if pointingField == 'NPF': PF_gal = [65.0, 30.0]
    if pointingField == 'SPF': PF_gal = [253.0, -30.0]

    PF_gal  = SkyCoord(PF_gal[0], PF_gal[1], frame='galactic', unit='deg')  # [deg]
    PF_icrs = PF_gal.icrs  # [deg]
    
    # Load PIC stars for each N-CAM visibility
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
    
    fig = plt.figure(figsize=(9,9))
    ax = plt.axes(projection='astro zoom', center=PF_icrs, radius='30 deg', rotate='184 deg')

    # Plot PIC1.1.0 stars after N-CAM visibility
    ax.plot(starPF06.ra.deg, starPF06.dec.deg, '.', c='skyblue',
            transform=ax.get_transform('world'), markersize=1, zorder=1)
    ax.plot(starPF12.ra.deg, starPF12.dec.deg, '.', c='deepskyblue',
            transform=ax.get_transform('world'), markersize=1, zorder=2)
    ax.plot(starPF18.ra.deg, starPF18.dec.deg, '.', c='dodgerblue',
            transform=ax.get_transform('world'), markersize=1, zorder=3)
    ax.plot(starPF24.ra.deg, starPF24.dec.deg, '.', c='royalblue',
            transform=ax.get_transform('world'), markersize=1, zorder=4)

    # Plot stars and add legend scaled to the stellar magnitudes
    if magStars is not None and len(magStars) > 0:
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
    #raGroups, decGroups = rf.getCameraGroupCoordinates(PF_icrs.ra.deg, PF_icrs.dec.deg, -8)
    #camPointing = SkyCoord(raGroups*u.deg, decGroups*u.deg, frame='icrs', unit='deg')  
    #ax.plot(camPointing.ra.deg, camPointing.dec.deg, 'r.', transform=ax.get_transform('world'), markersize=10, zorder=6)

    # Plot pointing of PIC1.1.0 and PIC2.0.0
    ax.plot(PF_icrs.ra.deg, PF_icrs.dec.deg, '*', transform=ax.get_transform('world'), ms=20, c='k', mfc='r', zorder=6)
    #ax.plot(277.18, 52.85, '*', transform=ax.get_transform('world'), ms=20, c='k', mfc='b', zorder=7)

    # Add-on's
    ax.scalebar((0.05, 0.05), 10 * u.deg).label()
    ax.compass(0.95, 0.05, 0.1)
    ax.grid(color='gray')

    # Settings
    if title is not None:
        ax.set_title(title, fontsize=fs+2)
    ax.set_xlabel('RA',  fontsize=fs)
    ax.set_ylabel('Dec', fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    ax.tick_params(axis='both', labelsize=fs)
    
    # Return figure
    return fig


#==============================================================#
#                         DISTRIBUTIONS                        #
#==============================================================#


def plotTeffvsRadius(fig0, starSample, title,
                     ds, df_dK, df_dG, df_dF,
                     sg, df_sgK, df_sgG, df_sgF,
                     df, ms_limit):

    """Distribution of Teff vs. Radius.

    Paramters
    ---------

    Return
    ------
    
    """
    
    # Fontsize
    ms = 3
    if starSample == 'P5': da = 0.2
    else: da = 0.0
    
    # Plots sub-giants
    plt.plot(df_sgK['Teff'], df_sgK['R'], 'o', alpha=0.5-da, markersize=ms, color='orange',      label=r'K$\,$IV')
    plt.plot(df_sgG['Teff'], df_sgG['R'], 'o', alpha=0.7-da, markersize=ms, color='greenyellow', label=r'G$\,$IV')
    plt.plot(df_sgF['Teff'], df_sgF['R'], 'o', alpha=0.7-da, markersize=ms, color='skyblue',     label=r'F$\,$IV')

    # Plot dwarfs
    plt.plot(df_dK['Teff'], df_dK['R'], 'o', alpha=0.4-da, markersize=ms, color='orangered', label=r'K$\,$V')
    plt.plot(df_dG['Teff'], df_dG['R'], 'o', alpha=0.4-da, markersize=ms, color='limegreen', label=r'G$\,$V')
    plt.plot(df_dF['Teff'], df_dF['R'], 'o', alpha=0.4-da, markersize=ms, color='royalblue', label=r'F$\,$V')

    # Plot selected targets
    plt.plot(df['Teff'], df['R'], 'o', mfc='none', mec='k', alpha=0.5, markersize=ms)

    # Compute main sequence devision
    dt = np.arange(np.min(ds['Teff']), np.max(ds['Teff']), 10)
    plt.plot(dt, ms_limit(dt), 'k-')

    # Settings
    plt.title(title)
    plt.xlabel(r'Effective temperature, $T_{\mathrm{eff}}$ [K]', fontsize=16)
    plt.ylabel(r'Stellar radius, $R$ [$R_{\odot}$]', fontsize=16)

    # Legend
    order = [3, 4, 5, 0, 1, 2]
    handles, labels = plt.gca().get_legend_handles_labels()
    h = [handles[idx] for idx in order]
    l = [labels[idx] for idx in order]
    plt.legend(h, l, ncol=2, loc='upper left', prop={'size':12},
               columnspacing=0.5, handletextpad=0)

    # Finito!
    return

    







    

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
    axes[0,0].set_xlabel(r'$P$ passband')
    #axes[0,0].set_xlabel(r'$V$ Johnson-Cousin')
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
    #axes[0,1].set_xlabel(r'$V$ Johnson-Cousin')
    axes[0,1].set_xlabel(r'$P$ passband')
    axes[0,1].set_ylabel('Number of stars')
    axes[0,1].tick_params(axis='x', which='minor', bottom=True, top=False)
    axes[0,1].tick_params(axis='x', which='major', bottom=True, top=False)
    axes[0,1].tick_params(axis='y', which='minor', left=False, right=False)
    axes[0,1].tick_params(axis='y', which='major', left=True, right=False)
    axes[0,1].grid(axis='y', color='gray', alpha=0.3)

    # Prepare bins and plot number distribution of contaminants per target

    numbinCon  = 1 + int(np.max(numConPerTar)/50)
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


#==============================================================#
#                     POINTING ERROR SOURCES                   #
#==============================================================#


def plotYawPitchRollTimeSeries(time, signals, units=["days", "arcsec"],
                               title=False, ylim=False, figsize=(10,10)):
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
        Axes matplotlib.pyplot handle object to be modified plot
    """

    # Datasets to loop over

    try:
        numData = len(signals)
    except:
        numData = 1

    # Handle yaxis limits

    if ylim is False:
        ylim = np.max(np.abs(signals))

    # Adjust linewidth after data

    if len(time) < 1e3: lw = 1
    else: lw = 0.5

    # Make plot

    fig = plt.figure(figsize=figsize)

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
        else:
            axes.set_xlabel('Time [{0}]'.format(units[0]))

        # Title

        if plot == 0: axes.set_title(title)

    # Adjust layout

    plt.tight_layout()
    plt.subplots_adjust(hspace = .001)

    # Finito!

    return axes





def plotYawPitchRollPSD(fig, time, signals, scale=1e-6, carbox=144,
                        title=False, labels=False, xmin=False, ylim=False, misreq=False):
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






def plotYawPitchRollJitter(time, signals, clabel, tpoint=100, lim=0.20,
                           cmap='gnuplot', plottype='short', title=False):

    """Pointing AOCS jitter illutrator. 

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
    None
    """

    # Hardcode valukes
    labels = ['Yaw [arcsec]', 'Pitch [arcsec]', 'Roll [arcsec]']
    nticks = 5
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
            im0 = ax[row, 0].scatter(signals[1][tpoint*row:tpoint*(row+1)],
                                     signals[0][tpoint*row:tpoint*(row+1)],
                                     c=time[tpoint*row:tpoint*(row+1)], s=sms, cmap=cmap, zorder=2)
            im1 = ax[row, 1].scatter(signals[2][tpoint*row:tpoint*(row+1)],
                                     signals[0][tpoint*row:tpoint*(row+1)],
                                     c=time[tpoint*row:tpoint*(row+1)], s=sms, cmap=cmap, zorder=2)
            im2 = ax[row, 2].scatter(signals[2][tpoint*row:tpoint*(row+1)],
                                     signals[1][tpoint*row:tpoint*(row+1)],
                                     c=time[tpoint*row:tpoint*(row+1)], s=sms, cmap=cmap, zorder=2)

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

        fig, ax = plt.subplots(1, 3, figsize=(10, 2.8))

        # Plot
        
        ax[0].plot(signals[1], signals[0], 'k-', alpha=al, lw=lw, zorder=1)
        ax[1].plot(signals[2], signals[0], 'k-', alpha=al, lw=lw, zorder=1)
        ax[2].plot(signals[2], signals[1], 'k-', alpha=al, lw=lw, zorder=1)
        im0 = ax[0].scatter(signals[1], signals[0], c=time, s=2, cmap='magma', zorder=2)
        im1 = ax[1].scatter(signals[2], signals[0], c=time, s=2, cmap='magma', zorder=2)
        im2 = ax[2].scatter(signals[2], signals[1], c=time, s=2, cmap='magma', zorder=2)

        # Labels
        
        ax[0].set_xlabel(labels[1])
        ax[0].set_ylabel(labels[0])
        ax[1].set_xlabel(labels[2])
        ax[1].set_ylabel(labels[0])
        ax[2].set_xlabel(labels[2])
        ax[2].set_ylabel(labels[1])

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

    # Finito!

    return fig, ax


#==============================================================#
#                        FREQUENCY ANALYSIS                    #
#==============================================================#


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


#==============================================================#
#                          PHOTOMETRY                          #
#==============================================================#


def plotPhotometry(df, time_unit=False, flux_unit=False, figsize=(8,5)):
    """
    PURPOSE: Function normalize the input flux and change time units to days. 

    NOTE: Function tailored to PLATOniums output format feather!

    PARAMETERS
    ----------
    data : pdarray
        Array containing at least a time and flux column (and potential flux_err)
    flux_err: boolen
       Wheather or not a flux_err column is present and should thus be added.

    RETURN
    ------
    data
    """

    # Create matplotlib object
    
    fig, ax = plt.subplots(1, 1, figsize=figsize)
        
    # Plot the input variable source
    if "flux_err" in df.columns:
        ax.errorbar(df["time"], df["flux"], yerr=df["flux_err"],
                    fmt=".", color='k', ecolor='darkgray', elinewidth=1,
                    capsize=0, alpha=aa, label="Raw flux", zorder=1)
    else:
        ax.plot(df["time"], df["flux"], 'k.', ms=ms, alpha=aa, label="Raw flux", zorder=1)

    # Plot a median filter
    
    if "flux_med" in df.columns:
        ax.plot(df["time"], df["flux_med"], 'b-', lw=lw, label='1h mdeian', zorder=2)
    
    # Show binned mean points if requested

    if "flux_bin" in df.columns:
        binsize = 1
        ax.plot(df["time"], df["flux_bin"], 'ro', ms=8, mec='k', label=f'{binsize}h bins', zorder=3)
            
    # Settings
    ax.set_xlim(df["time"].iloc[0], df["time"].iloc[-1])
    ax.set_xlabel(f"Time [{time_unit}]")
    ax.set_ylabel(f"Flux [{flux_unit}]")
    ax.legend(loc='best')
    
    return fig, ax



    


def plotNSRvsMagnitude(df, Ncam=1, tdur=3600., column=False, residuals=False,
                       camType="N", yscale="log", cmap="coolwarm", figsize=(10,6)):
    """
    PURPOSE:
    """
    
    # Create matplotlib object
    
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Set proper discrete cmap
    
    if column in ("group", "camera", "quarter", "ncam", "ncon", "flag"):
        # Fetch custom discrete colorbar used by matplotlib
        cbins = np.arange(df[column].min(), df[column].max()+2, 1)
        ticks = cbins + 0.5
        norm = discrete_colorbar(cbins=cbins, cmap=cmap)
    else:
        norm = None
    
    # Plot the input variable source
    
    if residuals:
        im = ax.scatter(df["mag"], df["res"], s=20, label="PlatoSim",
                        c=df[column], cmap=cmap, norm=norm)
        ax.set_ylabel('NSR Residuals [ppm]')
    elif column:
        im = ax.scatter(df["mag"], df["NSR"], s=20, label="PlatoSim",
                        c=df[column], cmap=cmap, norm=norm)
        ax.set_ylabel('NSR [ppm]')
    else:
        ax.plot(df["mag"], df["NSR"], 'k.', alpha=0.7, label="PlatoSim")
        ax.set_ylabel('NSR [ppm]')

    # Handle colorbar
    
    if norm is None:
        cb = plt.colorbar(im, extend="max", pad=0.01)
        cb.set_label(column)
    else:
        cb = plt.colorbar(im, extend="max", pad=0.01, spacing='proportional',
                          ticks=ticks, boundaries=cbins, format='%1i')
        cb.set_label(column)
        cb.minorticks_off()

    # Plot requirements
    
    if residuals == "camera":
        ax.axhline(y=108, c="darkorange", ls="--", label="AOCS camera req.: 108 ppm")
    elif residuals == "system":
        ax.axhline(y=9, c="red", ls="--", label="AOCS system req.: 9 ppm")
    elif residuals == "multi":
        ax.axhline(y=50, color="r", linestyle="-", label="24 N-CAM NSR req.: 50 ppm")
        ax.plot([11, 11, 11, 11], [88, 63, 52, 45], 'm*', ms=10, label="CBE of G0V 11 mag star")

    # Add grid

    # Calculate photon noise limit:
    # if not residuals:
    #     mag_range = np.linspace(df["mag"].min(), df["mag"].max(), 1000)
    #     NSR_cam01 = ut.getPhotonNoiseLimitNSR(mag_range, Ncam=1)
    #     NSR_cam06 = ut.getPhotonNoiseLimitNSR(mag_range, Ncam=6)
    #     NSR_cam12 = ut.getPhotonNoiseLimitNSR(mag_range, Ncam=12)
    #     NSR_cam18 = ut.getPhotonNoiseLimitNSR(mag_range, Ncam=18)
    #     NSR_cam24 = ut.getPhotonNoiseLimitNSR(mag_range, Ncam=24)
    #     NSR_phot = [NSR_cam01, NSR_cam06, NSR_cam12, NSR_cam18, NSR_cam24]
    #     cams = [1, 6, 12, 18, 24]
    #     for i, cam in zip(range(5), cams):
    #         ax.plot(mag_range, NSR_phot[i], '-', lw=2, label=f"{cam} N-CAM")

    # Force all yticks for log plot

    ax.set_yscale(yscale)
    subticks = [.1, .2, .3, .4, .5, .6, .7, .8, .9] 
    ax.yaxis.get_minor_locator().set_params(numticks=99, subs=subticks)
    
    # Settings
    
    ax.grid(color="lightgray")
    ax.set_xlabel(r'PLATO magnitude, $P$')
    ax.legend(loc='best')

    # Finito!
    
    return fig, ax

    
#==============================================================#
#                          ANIMATIONS                          #
#==============================================================#


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

    # Fetch StarPosition keys (i.e. "Exposure000000", etc)
    exposureGroupNames = list(f['StarPositions'].keys())[:-1] # -1 because Time is last column

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
            exposureGroupName = exposureGroupNames[imgNumber]
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
            plt.title(useTitle, fontsize=13)

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
        plt.xlabel('x [pixel]', fontsize=12)
        plt.ylabel('y [pixel]', fontsize=12)
            
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
