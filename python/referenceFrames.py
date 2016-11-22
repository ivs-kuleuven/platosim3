
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







def skyToFocalPlaneCoordinates(raStar, decStar, raOpticalAxis, decOpticalAxis, focalPlaneAngle, focalLength):

    """
    PURPOSE: Convert the equatorial sky coordinates (alpha, delta) of a star to focal plane coordinates (xFPprim, yFPprime),
             assuming a spherical pinhole camera.

    INPUT: raStar:           right ascension of the star [rad]
           decStar:          declination of the star [rad]
           raOpticalAxis:    right ascension of the optical axis [rad]
           decOpticalAxis:   declination of the optical axis [rad]
           focalPlaneAngle:  angle between the Y_FP axis and the Y'_FP axis: gamma_FP  [rad]
           focalLength:      focal length of the camera. Unit: see remarks.

    OUTPUT: xFPprime, yFPprime: cartesian coordinates of the project star in the focal plane in the FP-prime reference frame. 
                                Unit: see remarks:
    
    REMARK: - The unit of the cartesian focal plane coordinates is the same as the one of focalLength. 
              focalLength can be expressed in e.g. [mm] or in [pixels]. If focalLength == 1.0, then the corresponding
              focal plane coordinates are called "normalized coordinates."
    """

    # Convert the equatorial sky coordinate of the star to equatorial cartesian coordinates on the unit sphere

    xEQ = sin(pi/2.-decStar) * cos(raStar)
    yEQ = sin(pi/2.-decStar) * sin(raStar)
    zEQ = cos(pi/2.-decStar)

    vecEQ = array([xEQ, yEQ, zEQ])

    # Convert the equatorial cartesian coordinates to focal plane cartesian coordinates.

    rotMatrix1 = array([[ cos(raOpticalAxis), sin(raOpticalAxis), 0],      \
                        [-sin(raOpticalAxis), cos(raOpticalAxis), 0],      \
                        [                  0,                  0, 1]])

    rotMatrix2 = array([[sin(decOpticalAxis),  0, -cos(decOpticalAxis)],   \
                        [                  0,  1,                    0],   \
                        [cos(decOpticalAxis),  0,  sin(decOpticalAxis)]])

    vecFP = dot(rotMatrix2, dot(rotMatrix1, vecEQ))

    # Take into account the projection effect of the pinhole camera
    # Note that the pinhole reverses the image, hence the minus signs.

    xFP = -focalLength * vecFP[0]/vecFP[2]
    yFP = -focalLength * vecFP[1]/vecFP[2]

    # Convert the FP coordinates into FP' coordinates 

    xFPprime =  xFP * cos(focalPlaneAngle) + yFP * sin(focalPlaneAngle)
    yFPprime = -xFP * sin(focalPlaneAngle) + yFP * cos(focalPlaneAngle)

    # That's it!

    return xFPprime, yFPprime


def newSkyToFocalPlaneCoordinates(raStar, decStar, raPlatform, decPlatform, tiltAngle, azimuthAngle, focalPlaneAngle, focalLength):

    ## create the initial spacecraft coordinate system conversion matrix according 
    ## to PLATO-DLR-PL-TN-016_i1_2draft by Denis Griessbach

    ## the matrix transformations consist of rotations as well as translations, 
    ## which is why 4x4 matrices are used

    ## S/C reference frame

    ## the z-axis is the pointing coordinate of the spacecraft, the y-axis is lying in the equatorial plane 
    ## and the x-axis is calculated

    zSC = np.array([cos(decPlatform)*cos(raPlatform), cos(decPlatform)*sin(raPlatform), sin(decPlatform)])
    ySC = np.array([cos(raPlatform + pi/2), sin(raPlatform + pi/2), 0])
    xSC = np.cross(ySC, zSC)

    ## the first rotation matrix is the combination of the vectors

    rotSCtoEQ   = np.array([[xSC[0], ySC[0], zSC[0], 0.0], \
                            [xSC[1], ySC[1], zSC[1], 0.0], \
                            [xSC[2], ySC[2], zSC[2], 0.0], \
                            [  0.0 ,   0.0 ,   0.0 , 1.0]])

    ## include a rotation around the z - axis to align the x - axis to the average direction of the sun

    sunDirectionAngle = 0.0

    sinSun = sin(sunDirectionAngle)
    cosSun = cos(sunDirectionAngle)

    rotSun      = np.array([[cosSun, -sinSun, 0.0, 0.0], \
                            [sinSun,  cosSun, 0.0, 0.0], \
                            [  0.0 ,    0.0 , 1.0, 0.0], \
                            [  0.0 ,    0.0 , 0.0, 1.0]])
                       

    rotSCtoEQ = np.dot(rotSun, rotSCtoEQ)

    rotEQtoSC = np.transpose(rotSCtoEQ)


    ## the payload module reference frame

    ## the PLM frame coincides with the SC frame, except for a small mounting  offset, 
    ## which is a translation on the z-axis in positive direction

    distSCtoPLM = np.array([0.0, 0.0, 0.0])

    rotSCtoPLM  = np.array([[1.0, 0.0, 0.0, distSCtoPLM[0]], \
                            [0.0, 1.0, 0.0, distSCtoPLM[1]], \
                            [0.0, 0.0, 1.0, distSCtoPLM[2]], \
                            [0.0, 0.0, 0.0,     1.0       ]])


    ## the camera/FPA reference frame

    ## the x - axis should point towards the orientation screw, which can be achieved 
    ## by a rotation around the z  -axis

    screwAngle = 0.0

    sinScrew = sin(screwAngle)
    cosScrew = cos(screwAngle)

    rotScrew    = np.array([[cosScrew, -sinScrew, 0.0, 0.0], \
                            [sinScrew,  cosScrew, 0.0, 0.0], \
                            [   0.0  ,     0.0  , 1.0, 0.0], \
                            [   0.0  ,     0.0  , 0.0, 1.0]])

    ## the z - axis points towards the line of sight which is achieved by a 
    ## rotation of the z axis around the x - axis over the tilt - angle 

    sinTilt = sin(tiltAngle)
    cosTilt = cos(tiltAngle)

    rotTilt     = np.array([[1.0,   0.0  ,   0.0   , 0.0], \
                            [0.0, cosTilt, -sinTilt, 0.0], \
                            [0.0, sinTilt,  cosTilt, 0.0], \
                            [0.0,   0.0  ,   0.0   , 1.0]])



    ## the next step is the rotation over the focal plane angle

    sinOrient = sin(focalPlaneAngle)
    cosOrient = cos(focalPlaneAngle)

    rotOrient   = np.array([[cosOrient, -sinOrient, 0.0, 0.0], \
                            [sinOrient,  cosOrient, 0.0, 0.0], \
                            [   0.0   ,     0.0   , 1.0, 0.0], \
                            [   0.0   ,     0.0   , 0.0, 1.0]])


    rotPLMtoCAM = np.dot(rotOrient, np.dot(rotTilt, rotScrew))


    ## combine all transformation matrices

    rotEQtoCAM = np.dot(rotPLMtoCAM, np.dot(rotSCtoPLM, rotEQtoSC))

    ## create the cartesian coordinates of the star

    starEQ = np.array([cos(decStar)*cos(raStar),cos(decStar)*sin(raStar),sin(decStar), 1])

    ## calculation of the star coordinates in the CAM reference frame

    starCAM = np.dot(rotEQtoCAM, starEQ)

    ## calculate the mm values of the coordinates and return them

    xFPmm = - focalLength * starCAM[0]/starCAM[2]
    yFPmm = - focalLength * starCAM[1]/starCAM[2]

    return xFPmm, yFPmm



def focalPlaneToSkyCoordinates(xFPprime, yFPprime, raOpticalAxis, decOpticalAxis, focalPlaneAngle, focalLength):

    """
    PURPOSE: Convert the focal plane coordinates (xFPprim, yFPprime) of a star to equatorial sky coordinates (alpha, delta),
             assuming a spherical pinhole camera.

    INPUT: xFPprime:         Cartesian x-coordinate in the focal plane in the FP-prime reference frame [same unit as focalLength]
           yFPprime:         Cartesian y-coordinate in the focal plane in the FP-prime reference frame [same unit as focalLength]
           raOpticalAxis:    Right ascension of the optical axis [rad]
           decOpticalAxis:   Declination of the optical axis [rad]
           focalPlaneAngle:  Angle between the Y_FP axis and the Y'_FP axis: gamma_FP  [rad]
           focalLength:      focal length of the camera. Unit: [mm] or [pix] or [1.] or ...

    OUTPUT: raStar, decStar: Equatorial sky coordinates, right ascension and declination, of the star [rad]

    REMARK: The transformation assumes that the pinhole reverses the image.
    """

    # Convert the FP' coordinates in FP coordinates

    xFP =  xFPprime * cos(focalPlaneAngle) - yFPprime * sin(focalPlaneAngle)
    yFP =  xFPprime * sin(focalPlaneAngle) + yFPprime * cos(focalPlaneAngle)

    # Undo the reverse-image projection effect of the pinhole

    vecFP = array([-xFP/focalLength, -yFP/focalLength, 1.0])

    # Convert the focal plane cartesian coordinates to equatorial cartesian coordinates. 

    rotMatrix1 = array([[ sin(decOpticalAxis),  0, cos(decOpticalAxis)],    \
                        [                  0,   1,                    0],   \
                        [-cos(decOpticalAxis),  0, sin(decOpticalAxis)]])

    rotMatrix2 = array([[cos(raOpticalAxis), -sin(raOpticalAxis), 0],   \
                        [sin(raOpticalAxis),  cos(raOpticalAxis), 0],   \
                        [                  0,                  0, 1]])

    vecEQ = dot(rotMatrix2, dot(rotMatrix1, vecFP))


    # Convert the cartesian equatorial coordinates to equatorial sky coordinates

    norm = sqrt(vecEQ[0]*vecEQ[0] + vecEQ[1]*vecEQ[1] + vecEQ[2]*vecEQ[2]) 
    decStar = pi/2.0 - arccos(vecEQ[2]/norm);
    raStar = arctan2(vecEQ[1], vecEQ[0]);

    # Ensure that the right ascension is positive

    if isinstance(raStar, np.ndarray):
        raStar[raStar < 0.0] += 2.*pi
    else:
        if (raStar < 0.0):
            raStar += 2.*pi

    # That's it!

    return raStar, decStar












## \brief      Convert polar coordinates to cartesian coordinates
##
## \param[in]  distance  distance from the pole (reference point) [mm]
## \param[in]  angle     angle counter-clockwise from the x-axis [rad]
##
## \return     (xFPmm, yFPmm) Cartesian coordinates in the focal plane [mm]
##
def polarToCartesianFocalPlaneCoordinates(distance, angle):

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
def cartesianToPolarFocalPlaneCoordinates(xFPmm, yFPmm):
    
    angle = arctan2(yFPmm, xFPmm)      # [radians]
    distance = sqrt(xFPmm * xFPmm + yFPmm * yFPmm)

    return distance, angle









def undistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm):
    """
    PURPOSE:      Convert from undistorted to distorted focal plane coordinates
    
    INPUTS:       xFPmm  undistorted focal plane x-coordinate [mm]
                  yFPmm  undistorted focal plane y-coordinate [mm]
    
    OUTPUTS:      (xFPdist, yFPdist) distorted x and y coordinates [mm]
    """
    P = Polynomial(FIELD_DISTORTION["Coeff"])

    angle = arctan2(yFPmm, xFPmm)    # [radians]
    
    rFP = sqrt(xFPmm * xFPmm + yFPmm * yFPmm)
    rFPdist = P(rFP)
    
    xFPdist = cos(angle) * rFPdist
    yFPdist = sin(angle) * rFPdist
    
    return xFPdist, yFPdist










def distortedToUndistortedFocalPlaneCoordinates(xFPdist, yFPdist):
    """
    PURPOSE:     Convert from distorted to undistorted focal plane coordinates
    
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









def pixelToFocalPlaneCoordinates(xCCDpixel, yCCDpixel, pixelSize, ccdZeroPointX, ccdZeroPointY, CCDangle):

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











def focalPlaneToPixelCoordinates(xFPprime, yFPprime, pixelSize, ccdZeroPointX, ccdZeroPointY, CCDangle):

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

    cornersXmm, cornersYmm = pixelToFocalPlaneCoordinates(cornersXpix, cornersYpix, pixelSize, zeroPointXmm, zeroPointYmm, ccdAngle) 
    
    # That's it

    return cornersXmm, cornersYmm








def drawCCDsInSky(raOpticalAxis, decOpticalAxis, focalPlaneAngle, focalLength, pixelSize, normal=True):

    """
    PURPOSE: Project and plot the 4 CCDs of 1 camera on the sky

    INPUT: raOpticalAxis:   right ascension of the optical axis [rad]
           decOpticalAxis:  declination of the optical axis [rad]
           focalPlaneAngle: angle between the Y_FP axis and the Y'_FP axis: gamme_FP  [rad]
           focalLength:     [mm]
           pixelSize:       [micrometer]
           normal:         True for the normal camera configuration, False for the fast cameras

    OUTPUT: None

    TODO: Does not work yet for the fast cams
    """

    # Select the proper CCD codes depending on whether we're dealing with the nominal or the fast cams
    
    if normal == True:
        ccdCodes = ['A', 'B', 'C', 'D']
    else:
        ccdCodes = ['AF', 'BF', 'CF', 'DF']


    # Set up the colors to be used to draw each CCD. 
    # Different CCDs have different colors.

    color = {'A': 'b', 'AF': 'b', 'B': 'r', 'BF': 'r', 'C': 'g', 'CF': 'g', 'D': 'k', 'DF': 'k'}


    # Plot each of the 4 CCDs

    for ccdCode in ccdCodes:

        # Get the focal plane FP' coordinates of the CCD corners  [mm]

        cornersXmm, cornersYmm = computeCCDcornersInFocalPlane(ccdCode, pixelSize)

        # Compute the equatorial sky coordinates [rad] from the the focal plane FP' coordinates [mm] of the corners

        ra, dec = focalPlaneToSkyCoordinates(cornersXmm, cornersYmm, raOpticalAxis, decOpticalAxis, focalPlaneAngle, focalLength)

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











def drawCCDsInFocalPlane(pixelSize, normal=True):

    """
    PURPOSE: Plot the 4 CCDs in the focal plane in the FP' reference frame.
             May serve as a background to overplot the projected stars on the focal plane

    INPUT: pixelSize: size of 1 pixel [micron]
           normal: True for the normal camera configuration, False for the fast cameras

    OUTPUT: None

    """

    # Select the proper CCD codes depending on whether we're dealing with the nominal or the fast cams
    
    if normal == True:
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

    INPUT:    raStar:  right ascension of the star [rad]
              decStar: declination of the star [rad]
    
    OUTPUT:   Draw a red dot where the star is located on the CCD

    TODO: Update doc-string

    """

    normal = True  # FIXME: where can we specify that we use the fast or normal Camera

    if (sim["Camera/IncludeFieldDistortion"] == "yes")  or (sim["Camera/IncludeFieldDistortion"] == "1"):
        includeFieldDistortion = True
        FIELD_DISTORTION["Coeff"] = sim["Camera/FieldDistortion/Coefficients"]
        FIELD_DISTORTION["InverseCoeff"] = sim["Camera/FieldDistortion/InverseCoefficients"]
    else:
        includeFieldDistortion = False

    pixelSize = float(sim["CCD/PixelSize"])
    plateScale = float(sim["Camera/PlateScale"])
    raOpticalAxis = np.radians(float(sim["ObservingParameters/RApointing"]))
    decOpticalAxis = np.radians(float(sim["ObservingParameters/DecPointing"]))
    focalPlaneAngle = np.radians(float(sim["Camera/FocalPlaneOrientation"]))
    focalLength = float(sim["Camera/FocalLength"]) * 1000.0  # [m] -> [mm]
    ccdZeroPointX = float(sim["CCD/OriginOffsetX"])
    ccdZeroPointY = float(sim["CCD/OriginOffsetY"])
    ccdAngle = np.radians(float(sim["CCD/Orientation"]))

    xFPmm, yFPmm = skyToFocalPlaneCoordinates(raStar, decStar, raOpticalAxis, decOpticalAxis, focalPlaneAngle, focalLength)


    if includeFieldDistortion:
        xFPmm, yFPmm = undistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm)

    ccdCode, xCCD, yCCD = getCCDandPixelCoordinates(raStar, decStar, raOpticalAxis, decOpticalAxis, focalPlaneAngle, 
                                                    focalLength, plateScale, pixelSize, includeFieldDistortion, normal)

    if ccdCode == None:
        print ("Warning: DrawStarInFocalPlane(): The star doesn't fall on any of the CCDs.")
    else:
        drawPixelInFocalPlane(ccdCode, xCCD, yCCD, pixelSize)

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

    xFPprime, yFPprime = pixelToFocalPlaneCoordinates(xCCD, yCCD, pixelSize, zeroPointXmm, zeroPointYmm, ccdAngle)   # [mm]

    # Get the current axis

    currentAxis = plt.gca()
    currentAxis.add_patch( patches.Rectangle( (xFPprime, yFPprime), pixelSize / 1000.0, pixelSize / 1000.0, fill=False) )

    plt.plot(xFPprime, yFPprime, 'ro')
    plt.draw()

    return













def getCCDandPixelCoordinates(raStar, decStar, raPlatform, decPlatform, tiltAngle, azimuthAngle,  \
                              focalPlaneAngle, focalLength, plateScale, pixelSize, includeFieldDistortion, normal):

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
           normal:                 True for the normal camera configuration, False for the fast cameras

    OUTPUT: ccdCode: for normal camera: either 'A', 'B', 'C', or 'D'
                     for fast camer: either 'AF', 'BF', 'CF', 'DF'
                     if on no CCD: None    
            xCCDpix: x-coordinate (column number) of the star on the CCD  [pix]
                     if on no CCD: None
            yCCDpix: y-coordinate (row number) of the star on the CCD  [pix]
                     if on no CCD: None             
    """


    # Select the proper CCD codes depending on whether we're dealing with the nominal or the fast cams

    if normal == True:
        ccdCodes = ['A', 'B', 'C', 'D']
    else:
        ccdCodes = ['AF', 'BF', 'CF', 'DF']


    # Compute the (x,y) coordinates in the FP' reference system [mm]

    xFPmm, yFPmm = newSkyToFocalPlaneCoordinates(raStar, decStar, raPlatform, decPlatform, tiltAngle, azimuthAngle, focalPlaneAngle, focalLength)

    if includeFieldDistortion:
        xFPmm, yFPmm = undistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm)

    xCCDpixel = xFPmm / pixelSize * 1000.0
    yCCDpixel = yFPmm / pixelSize * 1000.0

    print(xCCDpixel, yCCDpixel)

    with open('/home/bert/PlatoSim3/inputfiles/test.txt', 'a') as fout:
      fout.write(str(float(xCCDpixel)) + "  " + str(float(yCCDpixel)) + "  " + str(float(xFPmm)) + "  " + str(float(yFPmm)) + '\n')


    # Find out if this falls on a CCD, and if yes which one.
    # Our approach: try each of the CCDs. Not elegant, but robust...  

    for ccdCode in ccdCodes:

        # Compute the position of the star in pixel coordinates, for the current CCD, 
        # disregarding the physical extend of the CCD

        zeroPointXmm = CCD[ccdCode]["zeroPointXmm"]
        zeroPointYmm = CCD[ccdCode]["zeroPointYmm"]
        ccdAngle     = CCD[ccdCode]["angle"]
        
        xCCDpix, yCCDpix = focalPlaneToPixelCoordinates(xFPmm, yFPmm, pixelSize, zeroPointXmm, zeroPointYmm, ccdAngle)

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













def calculateSubfieldAroundCoordinates(raStar, decStar, subfieldSizeX, subfieldSizeY, \
                                                             focalLength, plateScale, pixelSize, raPlatform, decPlatform, azimuthTelescope, tiltTelescope, 
                                                             focalPlaneAngle, includeFieldDistortion, normal):

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
             normal:                 True for the normal camera configuration, False for the fast cameras

    OUTPUTS: ccdCode: "A", "B", "C" or "D" if nominal=True, "AF", "BF", "CF" or "DF" otherwise
             xCCDpix: x-coordinate of the star in pixels (i.e. column number)
             yCCDpix: y-coordinate of the star in pixels (i.e. row number)

    REMARKS: - If the coordinates do not fall on any CCD, an error message is shown, followed by an exit(1)
             - If the star is too close to the edge for the given subfield size, and error message is shown,
               followed by an exit(1)
    """

    # Find out on which CCD the star falls, and the corresponding pixel coordinates

    ccdCode, xCCDpix, yCCDpix = getCCDandPixelCoordinates(raStar, decStar, raPlatform, decPlatform, tiltTelescope, azimuthTelescope, \
                                                          focalPlaneAngle, focalLength, plateScale, pixelSize, \
                                                          includeFieldDistortion, normal)

    # ccdCode, xCCDpix, yCCDpix = equatorialToFocalPlaneCoordinates(raStar, decStar, raPlatform, decPlatform, focalLength, azimuthTelescope, tiltTelescope, focalPlaneAngle, plateScale, pixelSize, includeFieldDistortion, normal)

    # If the CCD code is None, the star does not fall on any ccd -> error

    if ccdCode == None:
        print("Error: coordinates do not fall on any CCD")
        print("raStar, decStar = {0}, {1}".format(raStar, decStar))
        #print("Optical Axis (ra, dec, angle) = {0}, {1}, {2}". format(raOpticalAxis, decOpticalAxis, focalPlaneAngle))
        return None, None, None


    # If the star does fall on a CCD, check if it's not too close to the edge for the subfield to
    # be completely on the CCD.

    xCCDpix = round(xCCDpix)    # integer values
    yCCDpix = round(yCCDpix)
    firstRow = CCD[ccdCode]["firstRow"]     # different from nominal than for fast cams
    Ncols = CCD[ccdCode]["Ncols"]
    Nrows = CCD[ccdCode]["Nrows"]

    if     (xCCDpix - subfieldSizeX/2 < 0)        or (xCCDpix + subfieldSizeX/2 - 1 > Ncols-1)   \
        or (yCCDpix - subfieldSizeY/2 < firstRow) or (yCCDpix + subfieldSizeY/2 - 1> Nrows-1): 

        print("Error: pixel coordinates (row, col) = ({0},{1}) too close to the edge to accommodate subfield with size {2}x{3}" \
           .format(yCCDpix,xCCDpix,subfieldSizeX, subfieldSizeY))
        return None, None, None

    # That's it!

    return ccdCode, xCCDpix, yCCDpix











def setSubfieldAroundPixelCoordinates(sim, ccdCode, xCCDpixel, yCCDpixel, subfieldSizeX, subfieldSizeY):
    """
    PURPOSE:  Calculate the location of the subField such that it is centered on the star with 
              the given pixel coordinates.

    INPUTS:   sim            simulation for which the configuration file is adapted
              ccdCode:       for nominal camera: either 'A', 'B', 'C', 'D'
                             for fast camera: either 'AF', 'BF', 'CF', 'DF'
              xCCDpixel:     x-coordinate (column-number) of the star on the CCD  [pixel/float]
              yCCDpixel:     y-coordinate (row-number) of the star on the CCD  [pixel/float]
              subfieldSizeX: width (i.e. number of columns) of the sub-field [pixels]
              subfieldSizeY: height (i.e. number of rows) of the sub-field [pixels]

    OUTPUTS:  None
    """

    raStar, decStar = pixelToSkyCoordinates(sim, ccdCode, xCCDpixel, yCCDpixel)

    # TODO: determine nominal from the given ccdCode

    normal = True

    success = setSubfieldAroundCoordinates(sim, raStar, decStar, subfieldSizeX, subfieldSizeY, normal)

    if not success:
        print ("Warning: setSubfieldAroundPixelCoordinates() failed to set subField around the star.")

    return











def setSubfieldAroundCoordinates(sim, raStar, decStar, subfieldSizeX, subfieldSizeY, normal):
    
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
             normal:                 True for the normal camera configuration, False for the fast cameras

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
    focalPlaneAngle = np.deg2rad(float(sim["Camera/FocalPlaneOrientation"]))
    pixelSize = float(sim["CCD/PixelSize"]) 

    if (sim["Camera/IncludeFieldDistortion"] == "yes")  or (sim["Camera/IncludeFieldDistortion"] == "1"):
        includeFieldDistortion = True
    else:
        includeFieldDistortion = False


    # Derive the (RA, Dec) of the optical axis, given the (RA, Dec) of the platform, and the orientation
    # (azimith, tilt) of the telescope on the platform.

    #raOpticalAxis, decOpticalAxis = platformToTelescopePointingCoordinates(raPlatform, decPlatform, azimuthTelescope, tiltTelescope)    

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
                                                             focalLength, plateScale, pixelSize, raPlatform, decPlatform, azimuthTelescope, tiltTelescope, 
                                                             focalPlaneAngle, includeFieldDistortion, normal)
    
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

    if normal:
        sim["ObservingParameters/ExposureTime"] = 23
        sim["CCD/ReadoutTime"] = 2.5
    else:
        sim["ObservingParameters/ExposureTime"] = 2.3
        sim["CCD/ReadoutTime"] = 0.2
    
    # That's it

    return True












def skyToPixelCoordinates(sim, raStar, decStar, normal):
    """
    PURPOSE: Convert sky coordinates to pixel coordinates
    
    NOTE:   It is assumed that the configuration parameters in the sim object contains
            a correct (ra, dec) of the platform, a correct (azimuth, tilt) of the telescope,
            a valid value for the focal length, the plate scale, the pixel size, and that
            the switch to include distortion or not is set correctly.
            
    INPUT:  sim:        simulation for which the configuration file is adapted
            raStar:                 right ascension of the star [radians]
            decStar:                declination [radians]
            subfieldSizeX:          width (i.e. number of columns) of the subiield [pixels]
            subfieldSizeY:          height (i.e. number of rows) of the sub-field [pixels]
            normal:                 True for the normal camera configuration, False for the fast cameras
    
    OUTPUT: ccdCode
            xCCDpixel: column pixel coordinate of the star (real-valued) 
            yCCDpixel: row pixel coordinate of the star (real-valued)
    """
    
    if (sim["Camera/IncludeFieldDistortion"] == "yes")  or (sim["Camera/IncludeFieldDistortion"] == "1"):
        includeFieldDistortion = True
        FIELD_DISTORTION["Coeff"] = sim["Camera/FieldDistortion/Coefficients"]
        FIELD_DISTORTION["InverseCoeff"] = sim["Camera/FieldDistortion/InverseCoefficients"]
    else:
        includeFieldDistortion = False

    pixelSize = float(sim["CCD/PixelSize"])
    plateScale = float(sim["Camera/PlateScale"])
    focalLength = float(sim["Camera/FocalLength"]) * 1000.0       # [m] -> [mm]
    raPlatform = np.deg2rad(float(sim["ObservingParameters/RApointing"]))
    decPlatform = np.deg2rad(float(sim["ObservingParameters/DecPointing"]))
    focalPlaneAngle = float(sim["Camera/FocalPlaneOrientation"])
    azimuthTelescope = np.deg2rad(float(sim["Telescope/AzimuthAngle"]))
    tiltTelescope = np.deg2rad(float(sim["Telescope/TiltAngle"]))
    
    # Derive the (RA, Dec) of the optical axis, given the (RA, Dec) of the platform, and the orientation
    # (azimith, tilt) of the telescope on the platform.

    raOpticalAxis, decOpticalAxis = platformToTelescopePointingCoordinates(raPlatform, decPlatform, azimuthTelescope, tiltTelescope)    
    
    ccdCode, xCCDpixel, yCCDpixel = getCCDandPixelCoordinates(raStar, decStar, raOpticalAxis, decOpticalAxis, focalPlaneAngle, 
                                                              focalLength, plateScale, pixelSize, includeFieldDistortion, normal)

    return ccdCode, xCCDpixel, yCCDpixel
















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

    xFPmm, yFPmm = pixelToFocalPlaneCoordinates(xCCDpixel, yCCDpixel, pixelSize, ccdZeroPointX, ccdZeroPointY, ccdAngle)
    
    if includeFieldDistortion:
        xFPmm, yFPmm = distortedToUndistortedFocalPlaneCoordinates(xFPmm, yFPmm)
    
    ra, dec = focalPlaneToSkyCoordinates(xFPmm, yFPmm, raOpticalAxis, decOpticalAxis, focalPlaneAngle, focalLength)
    
    return ra, dec



## these are functions for an alternate way to calculate the position of the star or the spacecraft on the focal plane

def equatorialToFocalPlaneCoordinates(RA, dec, raPlatform, decPlatform, focalLength, azimuthTelescope, tiltTelescope, focalPlaneAngle, plateScale, pixelSize, includeFieldDistortion, normal):

    ## returns the values for the focal plane coordinates from the given right ascension and declination

    raOpticalAxis, decOpticalAxis = platformToTelescopePointingCoordinates(raPlatform, decPlatform, azimuthTelescope, tiltTelescope)

    ## construction of the z-axis

    zFP = np.array([cos(decPlatform)*cos(raPlatform), cos(decPlatform)*sin(raPlatform), sin(decPlatform)])

    ## the y axis lies in the equatorial plane

    yFP = np.array([cos(raPlatform + pi/2), sin(raPlatform + pi/2), 0])

    ## the x -axis is calculated as the cross product of the other two vectors

    xFP = np.cross(zFP, yFP)

    ## construction of the rotation matrix

    rotFpToEq = np.array([[xFP[0], yFP[0], zFP[0]], \
                          [xFP[1], yFP[1], zFP[1]], \
                          [xFP[2], yFP[2], zFP[2]]])

    rotEqToFp = np.transpose(rotFpToEq)

    ## convert the input parameters in a cartesian coordinate

    cartesianCoordinate = np.array([[cos(dec)*cos(RA)],[cos(dec)*sin(RA)],[sin(dec)]])

    ## rotation of the cartesian coordiante to the focal plane area

    focalPlaneCoordinate = np.dot(rotEqToFp, cartesianCoordinate)

    xFP = -focalLength * focalPlaneCoordinate[0]/focalPlaneCoordinate[2]
    yFP = -focalLength * focalPlaneCoordinate[1]/focalPlaneCoordinate[2]

    # Convert the FP coordinates into FP' coordinates 

    xFPmm =  xFP * cos(focalPlaneAngle) + yFP * sin(focalPlaneAngle)
    yFPmm = -xFP * sin(focalPlaneAngle) + yFP * cos(focalPlaneAngle)

    print (focalLength, focalPlaneAngle)

    if normal == True:
        ccdCodes = ['A', 'B', 'C', 'D']
    else:
        ccdCodes = ['AF', 'BF', 'CF', 'DF']


    if includeFieldDistortion:
        xFPmm, yFPmm = undistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm)

    # Find out if this falls on a CCD, and if yes which one.
    # Our approach: try each of the CCDs. Not elegant, but robust...

    xCCDpixel = xFPmm / pixelSize * 1000.0
    yCCDpixel = yFPmm / pixelSize * 1000.0

    with open('/home/bert/PlatoSim3/inputfiles/test.txt', 'a') as fout:
      fout.write(str(float(xCCDpixel)) + "  " + str(float(yCCDpixel)) + "  " + str(float(xFPmm)) + "  " + str(float(yFPmm)) + '\n')

    for ccdCode in ccdCodes:

        # Compute the position of the star in pixel coordinates, for the current CCD, 
        # disregarding the physical extend of the CCD

        zeroPointXmm = CCD[ccdCode]["zeroPointXmm"]
        zeroPointYmm = CCD[ccdCode]["zeroPointYmm"]
        ccdAngle     = CCD[ccdCode]["angle"]
        
        xCCDpix, yCCDpix = focalPlaneToPixelCoordinates(xFPmm, yFPmm, pixelSize, zeroPointXmm, zeroPointYmm, ccdAngle)

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


# def focalPlaneToEquatorialCoordinates(x, y, z):

#     ## returns the values for the focal plane coordinates from the given right ascension and declination

#     raPlatform = np.deg2rad(float(sim["ObservingParameters/RApointing"]))
#     decPlatform = np.deg2rad(float(sim["ObservingParameters/DecPointing"]))
#     azimuthTelescope = np.deg2rad(float(sim["Telescope/AzimuthAngle"]))
#     tiltTelescope = np.deg2rad(float(sim["Telescope/TiltAngle"]))

#     raOpticalAxis, decOpticalAxis = platformToTelescopePointingCoordinates(raPlatform, decPlatform, azimuthTelescope, tiltTelescope)

#     ## construction of the z-axis

#     zFP = np.array([cos(decOpticalAxis)*cos(raOpticalAxis), cos(decOpticalAxis)*sin(raOpticalAxis), sin(decOpticalAxis)])

#     ## the y axis lies in the equatorial plane

#     yFP = np.array([cos(raOpticalAxis + pi/2), sin(raOpticalAxis + pi/2), 0])

#     ## the x -axis is calculated as the cross product of the other two vectors

#     xFP = np.cross(zFP, yFP)

#     ## construction of the rotation matrix

#     rotFpToEq = np.array([[xFP[0], yFP[0], zFP[0]], \
#                            [xFP[1], yFP[1], zFP[1]], \
#                            [xFP[2], yFP[2], zFP[2]]])

#     ## rotation of the given coordinate to the equatorial coordinate system

#     equatorialCoordinate = dot(rotFpToEq, np.array([x],[y],[z])

#     ## calculating the right ascension and declination

#     norm = sqrt(equatorialCoordinate[0]*equatorialCoordinate[0]+equatorialCoordinate[1]*equatorialCoordinate[1]+equatorialCoordinate[2]*equatorialCoordinate[2])

#     delta = pi/2.0 - arccos(equatorialCoordinate[2]/norm)
#     alpha = arctan2(equatorialCoordinate[1], equatorialCoordinate[0])
#     if (alpha < 0.0): alpha += 2.0 * pi

#     return alpha, delta