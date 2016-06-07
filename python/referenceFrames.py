
import math
import numpy as np

from numpy import *
from numpy.polynomial import Polynomial

import matplotlib.pyplot as plt
import matplotlib.patches as patches

from matplotlib.path import Path





# CCD configuration
#
# ccdCode:      A,B,C,D: nominal cameras; AF, BF, CF, DF: fast cameras
# Ncols:        Number of exposed columns (column number varies along x-coordinate)
# Nrows:        Number of exposed rows (row number varies along y-coordinate)
# firstRow:     First row that is exposed. For the nominal cams this is simply row 0.
#               For the fast cams, the exposed rows are rows 2255 until 4510, which after the exposure 
#               are then frame-transfered to rows 0 until 2254.
# zeroPointXmm: x-coordinate of the (0,0) pixel of the CCD, in the FP' reference frame [mm]
# zeroPointYmm: y-coordinate of the (0,0) pixel of the CCD, in the FP' reference frame [mm]
# angle:        gamma_{ccd} [rad]: orientation angle of the CCD in the FP' reference frame

CCD = \
{
    'A'  : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 0,    'zeroPointXmm':  -1.0,   'zeroPointYmm': +82.162, 'angle': pi},
    'B'  : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 0,    'zeroPointXmm': +82.162, 'zeroPointYmm':  +1.0,   'angle': pi/2},
    'C'  : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 0,    'zeroPointXmm': -82.162, 'zeroPointYmm':  -1.0,   'angle': 3*pi/2},
    'D'  : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 0,    'zeroPointXmm':  +1.0,   'zeroPointYmm': -82.162, 'angle': 0},
    'AF' : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 2255, 'zeroPointXmm':  -1.0,   'zeroPointYmm': +82.162, 'angle': pi},
    'BF' : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 2255, 'zeroPointXmm': +82.162, 'zeroPointYmm':  +1.0,   'angle': pi/2},
    'CF' : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 2255, 'zeroPointXmm': -82.162, 'zeroPointYmm':  -1.0,   'angle': 3*pi/2},
    'DF' : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 2255, 'zeroPointXmm':  +1.0,   'zeroPointYmm': -82.162, 'angle': 0}
}

# FIELD DISTORTION
# 
# The field distortion is modelled by a 1D polynomial function of the 3rd degree.
# The FIELD_DISTORTION variable holds the coefficients of the polynomial and its inverse. The actual values
# are initialized by the setSubfieldAroundCoordinates() function when field distortion is requested.
# The evaluation is done by the Polynomial class of numpy.
# 
# Coeff:        Coefficients of the Polynomial
# InverseCoeff: Coefficients of the Inverse Polynomial
# 

FIELD_DISTORTION = {'Coeff' : [1.0, 0.0, 0.0, 0.0], 'InverseCoeff' : [-1.0, 0.0, 0.0, 0.0]}



def skyToAngularFocalPlaneCoordinates(raStar, decStar, raOpticalAxis, decOpticalAxis, focalPlaneAngle):

    """
    PURPOSE: computes the (x,y) coordinates in the focal plane of a star with given equatorial coordinates

    INPUT: raStar:          right ascension of the star [rad]
           decStar:         declination of the star [rad]
           raOpticalAxis:   right ascension of the optical axis [rad]
           decOpticalAxis:  declination of the optical axis [rad]
           focalPlaneAngle: angle between the Y_FP axis and the Y'_FP axis: gamme_FP  [rad]

    OUTPUT: xFPrad, yFPrad: Cartesian coordinate of the projected star in the focal plane in the FP-prime system [radians]
    """

    # Project the sky to the focal plane in the "FP" coordinate system

    denominator = cos(decOpticalAxis) * cos(decStar) * cos(raStar - raOpticalAxis) + sin(decOpticalAxis) * sin(decStar)
    xFP = (sin(decOpticalAxis) * cos(decStar) * cos(raStar - raOpticalAxis) - cos(decOpticalAxis) * sin(decStar)) / denominator
    yFP =  cos(decStar) * sin(raStar - raOpticalAxis) / denominator

    # Convert the FP coordinates into FP' coordinates 

    xFPrad =  xFP * cos(focalPlaneAngle) + yFP * sin(focalPlaneAngle)
    yFPrad = -xFP * sin(focalPlaneAngle) + yFP * cos(focalPlaneAngle)

    # Return the scaled coordinates

    return xFPrad, yFPrad







def angularFocalPlaneToSkyCoordinates(xFPrad, yFPrad, raOpticalAxis, decOpticalAxis, focalPlaneAngle):

    """
    PURPOSE: computes the (x,y) coordinates in the focal plane of a star with given equatorial coordinates

    INPUT: xFPrad            Cartesian coordinate of the projected star in the focal plane in the FP-prime system [radians]
           yFPrad            Cartesian coordinate of the projected star in the focal plane in the FP-prime system [radians]
           raOpticalAxis:    right ascension of the optical axis [rad]
           decOpticalAxis:   declination of the optical axis [rad]
           focalPlaneAngle:  angle between the Y_FP axis and the Y'_FP axis: gamme_FP  [rad]

    OUTPUT: raStar, decStar: Equatorial coordinates (right ascension and declination) of the star [rad]
    """

    xFP =  xFPrad * cos(focalPlaneAngle) - yFPrad * sin(focalPlaneAngle)
    yFP =  xFPrad * sin(focalPlaneAngle) + yFPrad * cos(focalPlaneAngle)

    if isscalar(xFP) and isscalar(yFP) and xFP == 0.0 and yFP == 0.0:
        return raOpticalAxis, decOpticalAxis

    rho = sqrt(xFP*xFP+yFP*yFP)
    c = arctan(rho)
    ra = raOpticalAxis + arctan2(yFP * sin(c), rho * cos(decOpticalAxis) * cos(c) + xFP * sin(decOpticalAxis) * sin(c))
    dec = arcsin(cos(c) * sin(decOpticalAxis) - (xFP * sin(c) * cos(decOpticalAxis)) / rho)

    return ra, dec








def angularToPlanarFocalPlaneCoordinates(xFPrad, yFPrad, focalLength):
    """
    PURPOSE: Convert from angular to planar focal plane coordinates, assuming no optical distortion.

    INPUT:   xFPrad:      Angular focal plane x-coordinate [rad]
             yFPrad:      Angular focal plane y-coordinate [rad]
             focalLength: Focal length of the telescope [mm]
    
    OUTPUT:  (xFPmm, yFPmm)    Planar focal plane x and y coordinates [mm]

    """

    xFPmm = tan(xFPrad) * focalLength
    yFPmm = tan(yFPrad) * focalLength

    return xFPmm, yFPmm








def planarToAngularFocalPlaneCoordinates(xFPmm, yFPmm, focalLength):
    """
    PURPOSE: Convert from planar to angulare focal plane coordinates, assuming no optical distortion.

    INPUT:   xFPmm:        Planar focal plane x-coordinate [mm]
             yFPmm:        Planar focal plane y-coordinate [mm]
             focalLength:  Focal length of the telescope [mm]
    
    OUTPUT:  (xFPrad, yFPrad)  Angular focal plane x and y coordinates [rad]

    """

    xFPrad = arctan(xFPmm / focalLength)
    yFPrad = arctan(yFPmm / focalLength)

    return xFPrad, yFPrad









## \brief      Convert polar coordinates to cartesian coordinates
##
## \param[in]  distance  distance from the pole (reference point) [mm]
## \param[in]  angle     angle counter-clockwise from the x-axis [rad]
##
## \return     (xFPmm, yFPmm) Cartesian coordinates in the focal plane [mm]
##
def polarToPlanarFocalPlaneCoordinates(distance, angle):

    xFPmm = cos(angle) * distance
    yFPmm = sin(angle) * distance

    return xFPmm, yFPmm







## \brief      Convert cartesian coordinates to polar coordinates
##
## \param[in]  xFPmm  x-axis cartesian coordinate in the focal plane [mm]
## \param[in]  yFPmm  y-axis cartesian coordinate in the focal plane [mm]
##
## \return     (distance, angle) polar coordinates in the focal plane
##
def planarToPolarFocalPlaneCoordinates(xFPmm, yFPmm):
    
    angle = arctan2(yFPmm, xFPmm)      # [radians]
    distance = sqrt(xFPmm * xFPmm + yFPmm * yFPmm)

    return distance, angle









def planarToDistortedFocalPlaneCoordinates(xFPmm, yFPmm):
    """
    PURPOSE:      Convert from planar to distorted focal plane coordinates
    
    INPUTS:       xFPmm  Planar focal plane x-coordinate [mm]
                  yFPmm  Planar focal plane y-coordinate [mm]
    
    OUTPUTS:      (xFPdist, yFPdist) distorted x and y coordinates [mm]
    """
    P = Polynomial(FIELD_DISTORTION["Coeff"])

    angle = arctan2(yFPmm, xFPmm)    # [radians]
    
    rFP = sqrt(xFPmm * xFPmm + yFPmm * yFPmm)
    rFPdist = P(rFP)
    
    xFPdist = cos(angle) * rFPdist
    yFPdist = sin(angle) * rFPdist
    
    return xFPdist, yFPdist










def distortedToPlanarFocalPlaneCoordinates(xFPdist, yFPdist):
    """
    PURPOSE:     Convert from distorted to planar focal plane coordinates
    
    INPUTS:      xFPdist  Distorted focal plane x-coordinate [mm]
                 yFPdist  DIstorted focal plane y-coordinate [mm]
    
    OUTPUTS:     (xFPmm, yFPmm) distorted x and y coordinates [mm]
    """
    IP = Polynomial(FIELD_DISTORTION["InverseCoeff"])
    
    angle = arctan2(yFPdist, xFPdist)  # [radians]
    
    rFP   = sqrt(xFPdist * xFPdist + yFPdist * yFPdist)
    rFPmm = IP(rFP)
    
    xFPmm = cos(angle) * rFPmm
    yFPmm = sin(angle) * rFPmm
    
    return xFPmm, yFPmm









def pixelToPlanarFocalPlaneCoordinates(xCCDpixel, yCCDpixel, pixelSize, ccdZeroPointX, ccdZeroPointY, CCDangle):

    """
    PUROSE: Given the (real-valued) pixel coordinates of the star on the CCD, compute the (x,y)
            coordinates in the FP' reference system (not the FP system!).

    INPUT: xCCDpixel     : x-coordinate (column-number) of the star on the CCD  [pixel]
           yCCDpixel     : y-coordinate (row-number) of the star on the CCD  [pixel]
           pixelSize     : size of 1 pixel in micron (not [mm]!)
           ccdZeroPointX : x-coordinate of the CCD (0,0) point in the FP' reference system [mm]
           ccdZeroPointY : y-coordinate of the CCD (0,0) point in the FP' reference system [mm]
           CCDangle      : CCD orientation angle in the FP' reference frame  [rad]

    OUTPUT: xFPprime: column pixel coordinate of the star (real-valued) [mm]
            yFPprime: row pixel coordinate of the star (real-valued) [mm]
    """

    # Convert the pixel coordinates into [mm] coordinates

    xCCDmm = xCCDpixel * pixelSize / 1000.0
    yCCDmm = yCCDpixel * pixelSize / 1000.0

    # Convert the CCD coordinates into FP' coordinates [mm]

    xFPprime = ccdZeroPointX + xCCDmm * cos(CCDangle) - yCCDmm * sin(CCDangle)
    yFPprime = ccdZeroPointY + xCCDmm * sin(CCDangle) + yCCDmm * cos(CCDangle)

    # That's it

    return xFPprime, yFPprime








def planarFocalPlaneToPixelCoordinates(xFPprime, yFPprime, pixelSize, ccdZeroPointX, ccdZeroPointY, CCDangle):

    """
    PUROSE: Compute the (real-valued) pixel coordinates of the star on the CCD, given the (x,y)
            coordinates in the FP' reference system (not the FP system!).

    INPUT: xFPprime      : x-coordinate of the star in the FP' reference system  [mm]
           yFPprime      : y-coordinate of the star in the FP' reference system  [mm]
           pixelSize     : size of 1 pixel [micron]  (not [mm]!)
           ccdZeroPointX : x-coordinate of the CCD (0,0) point in the FP' reference system [mm]
           ccdZeroPointY : y-coordinate of the CCD (0,0) point in the FP' reference system [mm]
           CCDangle      : CCD orientation angle in the FP' reference frame  [rad]

    OUTPUT: xCCDpixel: column pixel coordinate of the star (real-valued) 
            yCCDpixel: row pixel coordinate of the star (real-valued)
    """

    # Convert the FP' coordinates into CCD coordinates [mm]

    xCCDmm =  (xFPprime-ccdZeroPointX) * cos(CCDangle) + (yFPprime-ccdZeroPointY) * sin(CCDangle)
    yCCDmm = -(xFPprime-ccdZeroPointX) * sin(CCDangle) + (yFPprime-ccdZeroPointY) * cos(CCDangle)

    # Convert the [mm] coordinates into pixel coordinates

    xCCDpixel = xCCDmm / pixelSize * 1000.0
    yCCDpixel = yCCDmm / pixelSize * 1000.0

    # That's it

    return xCCDpixel, yCCDpixel








def computeCCDcornersInFocalPlane(ccdCode, pixelSize):

    """
    PURPOSE: Get the (x,y) coordinates of each of the 4 corners of the exposed part of the CCD
             in the FP' reference system

    INPUT: ccdCode:   one of the following: 'A', 'B', 'C', 'D', 'AF', 'BF', 'CF', 'DF'
           pixelSize: size of 1 pixel (micron)

    OUTPUT: cornersXmm: x-coordinates of each of the corners in the FP' reference system [mm]
            cornersYmm: y-coordinates of each of the corners in the FP' reference system [mm]
    """

    # Get the pixel coordinates of the 4 corners of the exposed part of the CCD
    # Note that the x-direction corresponds to the CCD columns, and the y-direction to the CCD rows.

    Nrows = CCD[ccdCode]["Nrows"]
    Ncols = CCD[ccdCode]["Ncols"]
    firstRow = CCD[ccdCode]["firstRow"]

    cornersXpix = array([0.0, Ncols-1, Ncols-1, 0.0])
    cornersYpix = array([firstRow, firstRow, Nrows-1, Nrows-1])
    
    # Convert to the x,y coordinates in the FP' reference frame

    zeroPointXmm = CCD[ccdCode]["zeroPointXmm"]
    zeroPointYmm = CCD[ccdCode]["zeroPointYmm"]
    ccdAngle     = CCD[ccdCode]["angle"]

    cornersXmm, cornersYmm = pixelToPlanarFocalPlaneCoordinates(cornersXpix, cornersYpix, pixelSize, zeroPointXmm, zeroPointYmm, ccdAngle) 
    
    # That's it

    return cornersXmm, cornersYmm








def drawCCDsInSky(raOpticalAxis, decOpticalAxis, focalPlaneAngle, focalLength, pixelSize, nominal=True):

    """
    PURPOSE: Project and plot the 4 CCDs of 1 camera on the sky

    INPUT: raOpticalAxis:   right ascension of the optical axis [rad]
           decOpticalAxis:  declination of the optical axis [rad]
           focalPlaneAngle: angle between the Y_FP axis and the Y'_FP axis: gamme_FP  [rad]
           focalLength:     [mm]
           pixelSize:       [micrometer]
           nominal:         True for the nominal camera configuration, False for the fast cameras

    OUTPUT: None

    TODO: Does not work yet for the fast cams
    """

    # Select the proper CCD codes depending on whether we're dealing with the nominal or the fast cams
    
    if nominal == True:
        ccdCodes = ['A', 'B', 'C', 'D']
    else:
        ccdCodes = ['AF', 'BF', 'CF', 'DF']


    # Set up the colors to be used to draw each CCD. 
    # Different CCDs have different colors.

    color = {'A': 'b', 'AF': 'b', 'B': 'r', 'BF': 'r', 'C': 'g', 'CF': 'g', 'D': 'k', 'DF': 'k'}


    # Plot each of the 4 CCDs

    for ccdCode in ccdCodes:

        # Get the planar FP' coordinates of the CCD corners  [mm]

        cornersXmm, cornersYmm = computeCCDcornersInFocalPlane(ccdCode, pixelSize)

        # Convert the planar FP' coordinates to angular FP' coordinates [rad]

        cornersXrad, cornersYrad = planarToAngularFocalPlaneCoordinates(cornersXmm, cornersYmm, focalLength)

        # Compute the equatorial sky coordinates [rad] from the the angular FP' coordinates [rad] of the corners

        ra, dec = angularFocalPlaneToSkyCoordinates(cornersXrad, cornersYrad, raOpticalAxis, decOpticalAxis, focalPlaneAngle)

        # Repeat the coordinates of the 1st corner, to plot a nice closed loop
        # Convert from radians to degrees

        ra  = append(ra, ra[0]) * 180 / pi
        dec = append(dec, dec[0]) * 180 / pi
        
        plt.plot(ra, dec, c=color[ccdCode])

        # Overplot the row closest to the readout register with a thicker line

        plt.plot([ra[0], ra[1]], [dec[0], dec[1]], c=color[ccdCode], linewidth=3)


    plt.xlabel("RA [deg]")
    plt.ylabel("Dec [deg]")
    plt.draw()

    # That's it

    return








def drawCCDsInFocalPlane(pixelSize, nominal=True):

    """
    PURPOSE: Plot the 4 CCDs in the focal plane in the FP' reference frame.
             May serve as a background to overplot the projected stars on the focal plane

    INPUT: pixelSize: size of 1 pixel [micron]
           nominal: True for the nominal camera configuration, False for the fast cameras

    OUTPUT: None

    """

    # Select the proper CCD codes depending on whether we're dealing with the nominal or the fast cams
    
    if nominal == True:
        ccdCodes = ['A', 'B', 'C', 'D']
    else:
        ccdCodes = ['AF', 'BF', 'CF', 'DF']


    # Set up the colors to be used to draw each CCD. 
    # Different CCDs have different colors.

    color = {'A': 'b', 'AF': 'b', 'B': 'r', 'BF': 'r', 'C': 'g', 'CF': 'g', 'D': 'k', 'DF': 'k'}


    # Plot each of the 4 CCDs

    for ccdCode in ccdCodes:

        # Get the corner coordinates in the FP' plane

        cornersXmm, cornersYmm = computeCCDcornersInFocalPlane(ccdCode, pixelSize)

        # Repeat the coordinates of the 1st corner, to plot a nice closed loop

        x = append(cornersXmm, cornersXmm[0])
        y = append(cornersYmm, cornersYmm[0])
        
        plt.plot(x, y, c=color[ccdCode])

        # Overplot the row closest to the readout register with a thicker line

        plt.plot([x[0], x[1]], [y[0], y[1]], c=color[ccdCode], linewidth=3)


    plt.xlabel("x_FP [mm]")
    plt.ylabel("y_FP [mm]")
    plt.draw()

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
    # disregarding the physical extend of the CCD

    zeroPointXmm = CCD[ccdCode]["zeroPointXmm"]
    zeroPointYmm = CCD[ccdCode]["zeroPointYmm"]
    ccdAngle     = CCD[ccdCode]["angle"]

    xFPprime, yFPprime = pixelToPlanarFocalPlaneCoordinates(xCCD, yCCD, pixelSize, zeroPointXmm, zeroPointYmm, ccdAngle)
    xFPprimeLL, yFPprimeLL = pixelToPlanarFocalPlaneCoordinates(xCCD - subfieldSizeX/2, yCCD - subfieldSizeY/2, \
        pixelSize, zeroPointXmm, zeroPointYmm, ccdAngle)
    xFPprimeUR, yFPprimeUR = pixelToPlanarFocalPlaneCoordinates(xCCD + subfieldSizeX/2, yCCD + subfieldSizeY/2, \
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

    INPUT:    raStar:  right ascension of the star [rad]
              decStar: declination of the star [rad]
    
    OUTPUT:   Draw a red dot where the star is located on the CCD

    TODO: Update doc-string

    """

    if (sim["Camera/IncludeFieldDistortion"] == "yes")  or (sim["Camera/IncludeFieldDistortion"] == "1"):
        includeFieldDistortion = True
        FIELD_DISTORTION["Coeff"] = sim["Camera/FieldDistortion/Coefficients"]
        FIELD_DISTORTION["InverseCoeff"] = sim["Camera/FieldDistortion/InverseCoefficients"]
    else:
        includeFieldDistortion = False

    pixelSize = float(sim["CCD/PixelSize"])
    raOpticalAxis = np.radians(float(sim["ObservingParameters/RApointing"]))
    decOpticalAxis = np.radians(float(sim["ObservingParameters/DecPointing"]))
    focalPlaneAngle = np.radians(float(sim["Camera/FocalPlaneOrientation"]))
    focalLength = float(sim["Camera/FocalLength"]) * 1000.0  # [m] -> [mm]
    ccdZeroPointX = float(sim["CCD/OriginOffsetX"])
    ccdZeroPointY = float(sim["CCD/OriginOffsetY"])
    ccdAngle = np.radians(float(sim["CCD/Orientation"]))

    xFPrad, yFPrad = skyToAngularFocalPlaneCoordinates(raStar, decStar, raOpticalAxis, decOpticalAxis, focalPlaneAngle)
    xFPmm, yFPmm = angularToPlanarFocalPlaneCoordinates(xFPrad, yFPrad, focalLength)
    if includeFieldDistortion:
        xFPmm, yFPmm = planarToDistortedFocalPlaneCoordinates(xFPmm, yFPmm)

    xCCD, yCCD = planarFocalPlaneToPixelCoordinates(xFPmm, yFPmm, pixelSize, ccdZeroPointX, ccdZeroPointY, ccdAngle)

    # TODO: Determine the ccdCode

    drawPixelInFocalPlane("B", xCCD, yCCD, pixelSize)

    return









def drawPixelInFocalPlane(ccdCode, xCCD, yCCD, pixelSize):

    """
    PURPOSE: Plot a pixel from a particular CCD in the focal plane. The actual position in millimeter
             is shown as a red dot, while the pixel itself is drawn as a rectangle with edge pixelSize.

    INPUTS:  ccdCode:   for nominal camera: either 'A', 'B', 'C', or 'D'
                        for fast camer: either 'AF', 'BF', 'CF', 'DF'
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

    xFPprime, yFPprime = pixelToPlanarFocalPlaneCoordinates(xCCD, yCCD, pixelSize, zeroPointXmm, zeroPointYmm, ccdAngle)

    # Get the current axis

    currentAxis = plt.gca()
    currentAxis.add_patch( patches.Rectangle( (xFPprime, yFPprime), pixelSize / 1000.0, pixelSize / 1000.0, fill=False) )

    plt.plot(xFPprime, yFPprime, 'ro')
    plt.draw()

    return









def getCCDandPixelCoordinates(raStar, decStar, raOpticalAxis, decOpticalAxis, focalPlaneAngle,  \
                              focalLength, plateScale, pixelSize, includeFieldDistortion=True, nominal=True):

    """
    PURPOSE: Given the equatorial coordinates of a star, find out on which CCD it falls ('A', 'B', ...)
             and compute the pixel coordinates of the star on this CCD.

    INPUT: raStar:                 right ascension of the star [rad]
           decStar:                declination of the star [rad]      
           raOpticalAxis:          right ascension of the optical axis [rad]
           decOpticalAxis:         declination of the optical axis [rad]
           focalPlaneAngle:        angle between the Y_FP axis and the Y'_FP axis: gamme_FP  [rad]
           focalLength:            focal length of the telescope [mm]
           plateScale:             [arcsec/micron]
           pixelSize:              [micrometer]
           includeFieldDistortion: True to include field distortion in coordinate transformations, false otherwise
           nominal:                True for the nominal camera configuration, False for the fast cameras

    OUTPUT: ccdCode: for nominal camera: either 'A', 'B', 'C', or 'D'
                     for fast camer: either 'AF', 'BF', 'CF', 'DF'
                     if on no CCD: None    
            xCCDpix: x-coordinate (column number) of the star on the CCD  [pix]
                     if on no CCD: None
            yCCDpix: y-coordinate (row number) of the star on the CCD  [pix]
                     if on no CCD: None             
    """


    # Select the proper CCD codes depending on whether we're dealing with the nominal or the fast cams

    if nominal == True:
        ccdCodes = ['A', 'B', 'C', 'D']
    else:
        ccdCodes = ['AF', 'BF', 'CF', 'DF']


    # Compute the (x,y) coordinates in the FP' reference system [mm]

    xFPrad, yFPrad = skyToAngularFocalPlaneCoordinates(raStar, decStar, raOpticalAxis, decOpticalAxis, focalPlaneAngle)
    xFPmm, yFPmm = angularToPlanarFocalPlaneCoordinates(xFPrad, yFPrad, focalLength)

    if includeFieldDistortion:
        xFPmm, yFPmm = planarToDistortedFocalPlaneCoordinates(xFPmm, yFPmm)

    # Find out if this falls on a CCD, and if yes which one.
    # Our approach: try each of the CCDs. Not elegant, but robust...

    for ccdCode in ccdCodes:

        # Compute the position of the star in pixel coordinates, for the current CCD, 
        # disregarding the physical extend of the CCD

        zeroPointXmm = CCD[ccdCode]["zeroPointXmm"]
        zeroPointYmm = CCD[ccdCode]["zeroPointYmm"]
        ccdAngle     = CCD[ccdCode]["angle"]
        
        xCCDpix, yCCDpix = planarFocalPlaneToPixelCoordinates(xFPmm, yFPmm, pixelSize, zeroPointXmm, zeroPointYmm, ccdAngle)

        # Check if the star falls on the exposed area of the CCD. If not: go to next CCD

        Nrows = CCD[ccdCode]["Nrows"]
        Ncols = CCD[ccdCode]["Ncols"]
        firstRow = CCD[ccdCode]["firstRow"]

        if (xCCDpix < 0) or (yCCDpix < firstRow): continue
        if (xCCDpix >= Ncols) or (yCCDpix >= Nrows): continue

        # If we arrive here, we found a CCD on which the star is located

        return ccdCode, xCCDpix, yCCDpix


    # If we arrive here, the star does not fall on any CCD  

    return None, None, None









def platformToTelescopePointingCoordinates(alphaPlatform, deltaPlatform, azimuthAngle, tiltAngle):

    """
    PURPOSE: Given the platform pointing coordinates (i.e. the sky coordinates of the jitter axis)
             and the orientation of the telescope on the platform, compute the sky coordinates of
             the optical axis of the telescope.
             See also: PLATO-KUL-PL-TN-001

    INPUT: alphaPlatform:  Right Ascension of the Platform pointing axis [rad]
           deltaPlatform:  Declination of the platfor pointing axis [rad]
           azimuthAngle:   azimuth of the telescope on the platform [rad]
           tiltAngle:      tilt angle between platform and telescope pointing axes [rad]

    OUTPUT: alphaTelescope: right ascension of the optical axis of the telescope [rad]
            deltaTelescope: declination of the optical axis of the telescope [rad]

    """

    # Specify the coordinates in the equatorial reference frame of the unit vector zSC,
    # corresponding to the jitter (Z) axis of the SpaceCraft

    zSC = array([cos(alphaPlatform)*cos(deltaPlatform), sin(alphaPlatform)*cos(deltaPlatform), sin(deltaPlatform)])

    # Construct the rotation matrix for a rotation around the zSC axis over the azimuth angle
    # {ux, uy, uz} is short-hand notation.

    ux = zSC[0]
    uy = zSC[1]
    uz = zSC[2]

    cosAngle = cos(azimuthAngle)
    sinAngle = sin(azimuthAngle)

    rotAzimuth = array([[cosAngle+ux*ux*(1-cosAngle),    ux*uy*(1-cosAngle)-uz*sinAngle, ux*uz*(1-cosAngle)+uy*sinAngle], \
                        [uy*ux*(1-cosAngle)+uz*sinAngle, cosAngle+uy*uy*(1-cosAngle),    uy*uz*(1-cosAngle)-ux*sinAngle], \
                        [uz*ux*(1-cosAngle)-uy*sinAngle, uz*uy*(1-cosAngle)+ux*sinAngle, cosAngle+uz*uz*(1-cosAngle)]])


    # The goal of the rotZ rotation matrix is to rotate the ySC unit vector (corresponding to the y-axis in 
    # the spacecraft reference frame) in the azimuth direction of the telescope. Rather than using ySC, we use
    # another reference vector yRef as defined below. y0 is perpendicular to zSC, and has the advantage that the 
    # exact orientation of the spacecraft (i.e. in which direction the sunshield is pointing) is not needed.

    yRef = array([-sin(alphaPlatform), cos(alphaPlatform), 0.0])

    # Rotate this reference vector
    # Note: numpy.matrix uses algebraic multiplication.

    yAzimuth = dot(rotAzimuth, yRef)

    # Next, construct the rotation matrix for a rotation around the yAzimuth vector over the tilt angle
    # of the telescope. The tilt angle is the angle between the optical axis of the telescope and the
    # the jitter Z-axis of the platform.

    ux = yAzimuth[0]
    uy = yAzimuth[1]
    uz = yAzimuth[2]

    cosAngle = cos(tiltAngle)
    sinAngle = sin(tiltAngle)    

    rotTilt = array([[cosAngle+ux*ux*(1-cosAngle),    ux*uy*(1-cosAngle)-uz*sinAngle, ux*uz*(1-cosAngle)+uy*sinAngle], \
                     [uy*ux*(1-cosAngle)+uz*sinAngle, cosAngle+uy*uy*(1-cosAngle),    uy*uz*(1-cosAngle)-ux*sinAngle], \
                     [uz*ux*(1-cosAngle)-uy*sinAngle, uz*uy*(1-cosAngle)+ux*sinAngle, cosAngle+uz*uz*(1-cosAngle)]])


    # Compute the unit vector zOA in the direction of the telescope's optical axis

    zOA = dot(rotTilt, zSC);


    # zOA now contains the cartesian coordinates of the optical axis in the equatorial reference frame. 
    # Compute the equatorial sky coordinates [rad] from the cartesian coordinates.

    norm = sqrt(zOA[0]*zOA[0]+zOA[1]*zOA[1]+zOA[2]*zOA[2])

    deltaTelescope = pi/2.0 - arccos(zOA[2]/norm)
    alphaTelescope = arctan2(zOA[1], zOA[0])
    if (alphaTelescope < 0.0): alphaTelescope += 2.0 * pi

    return alphaTelescope, deltaTelescope













def calculateSubfieldAroundCoordinates(raStar, decStar, subfieldSizeX, subfieldSizeY, focalLength, plateScale, pixelSize, \
                                       raOpticalAxis, decOpticalAxis, focalPlaneAngle, includeFieldDistortion=True, nominal=True):

    """
    PURPOSE: Calculates the location of the subfield such that the star with coordinates (raStar, decStar)
             is centered in the subfield. This function is used by setSubfieldAroundCoordinates() and usually
             does not need to be called by the user.

    NOTE: This function requires (raOpticalAxis, decOpticalAxis) while the function setSubfieldAroundCoordinates()
          requires (raPlatform, decPlatform).

    INPUTS:  raStar:                 right ascension [rad]
             decStar:                declination [rad]
             subfieldSizeX:          full width (# of columns) of the subfield [pix]
             subfieldSizeY:          full height (#of rows) of the subfield [pix]
             focalLength:            focal length of the telescope [mm]
             plateScale:             [arcsec/micron]
             pixelSize:              [micrometer]
             raOpticalAxis:          right ascension of the optical axis [rad]
             decOpticalAxis:         declination of the optical axis [rad]
             focalPlaneAngle:        angle between the Y_FP axis and the Y'_FP axis: gamme_FP  [rad]
             includeFieldDistortion: True to include field distortion in coordinate transformations, false otherwise
             nominal:                True for the nominal camera configuration, False for the fast cameras

    OUTPUTS: ccdCode: "A", "B", "C" or "D" if nominal=True, "AF", "BF", "CF" or "DF" otherwise
             xCCDpix: x-coordinate of the star in pixels (i.e. column number)
             yCCDpix: y-coordinate of the star in pixels (i.e. row number)

    REMARKS: - If the coordinates do not fall on any CCD, an error message is shown, followed by an exit(1)
             - If the star is too close to the edge for the given subfield size, and error message is shown,
               followed by an exit(1)
    """

    # Find out on which CCD the star falls, and the corresponding pixel coordinates

    ccdCode, xCCDpix, yCCDpix = getCCDandPixelCoordinates(raStar, decStar, raOpticalAxis, decOpticalAxis, \
        focalPlaneAngle, focalLength, plateScale, pixelSize, includeFieldDistortion, nominal)

    # If the CCD code is None, the star does not fall on any ccd -> error

    if ccdCode == None:
        print("Error: coordinates do not fall on any CCD")
        print("raStar, decStar = {0}, {1}".format(raStar, decStar))
        print("Optical Axis (ra, dec, angle) = {0}, {1}, {2}". format(raOpticalAxis, decOpticalAxis, focalPlaneAngle))
        return None, None, None


    # If the star does fall on a CCD, check if it's not too close to the edge for the subfield to
    # be completely on the CCD.

    xCCDpix = round(xCCDpix)    # integer values
    yCCDpix = round(yCCDpix)
    firstRow = CCD[ccdCode]["firstRow"]     # different from nominal than for fast cams
    Ncols = CCD[ccdCode]["Ncols"]
    Nrows = CCD[ccdCode]["Nrows"]

    if     (xCCDpix - subfieldSizeX/2 < 0)        or (xCCDpix + subfieldSizeX/2 > Ncols-1)   \
        or (yCCDpix - subfieldSizeY/2 < firstRow) or (yCCDpix + subfieldSizeY/2 > Nrows-1): 

        print("Error: pixel coordinates (row, col) = ({0},{1}) too close to the edge to accommodate subfield with size {2}x{3}" \
           .format(yCCDpix,xCCDpix,subfieldSizeX, subfieldSizeY))
        return None, None, None

    # That's it!

    return ccdCode, xCCDpix, yCCDpix












def setSubfieldAroundCoordinates(sim, raStar, decStar, subfieldSizeX, subfieldSizeY, nominal=True):
    
    """
    PURPOSE: Calculates the location of the sub-field such that it is centred on the star 
             with the given sky coordinates.  Depending on the CCD (in nomincal mode:
             "A", "B", "C", or "D"; in fast mode: "AF", "BF", "CF", or "DF"), the 
             configuration file for the given simulation are adapted.  These include the
             pre-defined CCD position, the dimensions of the CCD (and also of the sub-field,
             although this is not affected by the calculations), the sub-field zeropoint
             and the exposure time. 

    NOTE: This function calls the calculateSubfieldAroundCoordinates() function.

    NOTE: It is assumed that the configuration parameters in the sim object contains
          a correct (ra, dec)  of the platform, a correct (azimuth, tilt) of the telescope,
          a valid values for the focal length, the plate scale, the pixel size, and that
          the switch to include distortion or not is set correctly

    INPUTS:  sim:                    simulation for which the configuration file is adapted
             raStar:                 right ascension of the star [radians]
             decStar:                declination [radians]
             subfieldSizeX:          width (i.e. number of columns) of the subiield [pixels]
             subfieldSizeY:          height (i.e. number of rows) of the sub-field [pixels]
             nominal:                True for the nominal camera configuration, False for the fast cameras

    OUTPUT: True if the CCD code (i.e. the pre-defined CCD position) could be determined, False otherwise 

    REMARKS: - If the coordinates do not fall on any CCD, an error message is shown, followed by an exit(1)
             - If the star is too close to the edge for the given subfield size, and error message is shown,
               followed by an exit(1)
    """
    

    # Find out some instrumental characteristics from the sim object

    raPlatform = np.deg2rad(float(sim["ObservingParameters/RApointing"]))
    decPlatform = np.deg2rad(float(sim["ObservingParameters/DecPointing"]))
    azimuthTelescope = np.deg2rad(float(sim["Telescope/AzimuthAngle"]))
    tiltTelescope = np.deg2rad(float(sim["Telescope/TiltAngle"]))
    focalLength = float(sim["Camera/FocalLength"]) * 1000.0       # [m] -> [mm]
    plateScale = float(sim["Camera/PlateScale"])          
    focalPlaneAngle = float(sim["Camera/FocalPlaneOrientation"])
    pixelSize = float(sim["CCD/PixelSize"]) 

    if (sim["Camera/IncludeFieldDistortion"] == "yes")  or (sim["Camera/IncludeFieldDistortion"] == "1"):
        includeFieldDistortion = True
    else:
        includeFieldDistortion = False


    # Derive the (RA, Dec) of the optical axis, given the (RA, Dec) of the platform, and the orientation
    # (azimith, tilt) of the telescope on the platform.

    raOpticalAxis, decOpticalAxis = platformToTelescopePointingCoordinates(raPlatform, decPlatform, azimuthTelescope, tiltTelescope)    

    # When the user requested to include field distortion, update the Simulation input parameter and
    # initialize the field distortion global that will be used by the distortion functions.

    if includeFieldDistortion:
        sim["Camera/IncludeFieldDistortion"] = "yes"

        FIELD_DISTORTION["Coeff"] = sim["Camera/FieldDistortion/Coefficients"]
        FIELD_DISTORTION["InverseCoeff"] = sim["Camera/FieldDistortion/InverseCoefficients"]
    else:
        sim["Camera/IncludeFieldDistortion"] = "no"


    # Compute the position of the subfield.
    # xPix and yPix are the CCD coordinates of the star, given a 4510x4510 CCD [colNumber, rowNumber].

    ccdCode, xPix, yPix = calculateSubfieldAroundCoordinates(raStar, decStar, subfieldSizeX, subfieldSizeY, \
            focalLength, plateScale, pixelSize, raOpticalAxis, decOpticalAxis, focalPlaneAngle, includeFieldDistortion, nominal)
    
    if ccdCode == None:
        return False
    
    CCDSizeX         = CCD[ccdCode]["Ncols"]
    CCDSizeY         = CCD[ccdCode]["Nrows"]
    CCDOriginOffsetX = CCD[ccdCode]["zeroPointXmm"]
    CCDOriginOffsetY = CCD[ccdCode]["zeroPointYmm"]
    CCDOrientation   = CCD[ccdCode]["angle"]

    # If we arrive here, there is no problem accommodating the entire sufield on the CCD

    sim["CCD/OriginOffsetX"] = str(CCDOriginOffsetX)
    sim["CCD/OriginOffsetY"] = str(CCDOriginOffsetY)
    sim["CCD/Orientation"] = str(rad2deg(CCDOrientation))

    sim["CCD/NumColumns"] = CCDSizeX
    sim["CCD/NumRows"] = CCDSizeY

    sim["SubField/ZeroPointRow"] = str(int(yPix - subfieldSizeY/2))
    sim["SubField/ZeroPointColumn"] = str(int(xPix - subfieldSizeX/2))
    sim["SubField/NumRows"] = str(subfieldSizeY)
    sim["SubField/NumColumns"] = str(subfieldSizeX)

    # Set the exposure and the readout time, depending on fast vs nominal cams

    if nominal:
        sim["ObservingParameters/ExposureTime"] = 23
        sim["CCD/ReadoutTime"] = 2.5
    else:
        sim["ObservingParameters/ExposureTime"] = 2.3
        sim["CCD/ReadoutTime"] = 0.2
    
    # That's it

    return True












def pixelToSkyCoordinates(sim, ccdCode, xCCDpixel, yCCDpixel):
    """
    PURPOSE: Convert pixel coordinates to sky coordinates
    
    NOTE:   It is assumed that the configuration parameters in the sim object contains
            a correct (ra, dec) of the platform, a correct (azimuth, tilt) of the telescope,
            a valid value for the focal length, the plate scale, the pixel size, and that
            the switch to include distortion or not is set correctly.
            
    INPUT:  sim:        simulation for which the configuration file is adapted
            ccdCode:    for nominal camera: either 'A', 'B', 'C', 'D'
                        for fast camera: either 'AF', 'BF', 'CF', 'DF'
            xCCDpixel:  x-coordinate (column-number) of the star on the CCD  [pixel]
            yCCDpixel:  y-coordinate (row-number) of the star on the CCD  [pixel]
    
    OUTPUT: raStar, decStar: Equatorial coordinates (right ascension and declination) of the star [rad]
    """
    
    if (sim["Camera/IncludeFieldDistortion"] == "yes")  or (sim["Camera/IncludeFieldDistortion"] == "1"):
        includeFieldDistortion = True
        FIELD_DISTORTION["Coeff"] = sim["Camera/FieldDistortion/Coefficients"]
        FIELD_DISTORTION["InverseCoeff"] = sim["Camera/FieldDistortion/InverseCoefficients"]
    else:
        includeFieldDistortion = False

    pixelSize = float(sim["CCD/PixelSize"])
    focalLength = float(sim["Camera/FocalLength"]) * 1000.0       # [m] -> [mm]
    raPlatform = np.deg2rad(float(sim["ObservingParameters/RApointing"]))
    decPlatform = np.deg2rad(float(sim["ObservingParameters/DecPointing"]))
    focalPlaneAngle = float(sim["Camera/FocalPlaneOrientation"])
    azimuthTelescope = np.deg2rad(float(sim["Telescope/AzimuthAngle"]))
    tiltTelescope = np.deg2rad(float(sim["Telescope/TiltAngle"]))

    ccdZeroPointX = CCD[ccdCode]['zeroPointXmm']
    ccdZeroPointY = CCD[ccdCode]['zeroPointYmm']
    ccdAngle = CCD[ccdCode]['angle']
    
    # Derive the (RA, Dec) of the optical axis, given the (RA, Dec) of the platform, and the orientation
    # (azimith, tilt) of the telescope on the platform.

    raOpticalAxis, decOpticalAxis = platformToTelescopePointingCoordinates(raPlatform, decPlatform, azimuthTelescope, tiltTelescope)    

    xFPmm, yFPmm = pixelToPlanarFocalPlaneCoordinates(xCCDpixel, yCCDpixel, pixelSize, ccdZeroPointX, ccdZeroPointY, ccdAngle)
    
    if includeFieldDistortion:
        xFPmm, yFPmm = distortedToPlanarFocalPlaneCoordinates(xFPmm, yFPmm)
    
    xFPrad, yFPrad = planarToAngularFocalPlaneCoordinates(xFPmm, yFPmm, focalLength)
    ra, dec = angularFocalPlaneToSkyCoordinates(xFPrad, yFPrad, raOpticalAxis, decOpticalAxis, focalPlaneAngle)
    
    return ra, dec

