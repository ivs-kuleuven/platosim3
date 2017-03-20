
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
    'A'  : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 0,    'zeroPointXmm':  -1.0, 'zeroPointYmm': +82.162, 'angle': pi},
    'B'  : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 0,    'zeroPointXmm':  -1.0, 'zeroPointYmm': +82.162, 'angle': 3*pi/2},
    'C'  : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 0,    'zeroPointXmm':  -1.0, 'zeroPointYmm': +82.162, 'angle': pi/2},
    'D'  : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 0,    'zeroPointXmm':  -1.0, 'zeroPointYmm': +82.162, 'angle': 0},
    'AF' : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 2255, 'zeroPointXmm':  -1.0, 'zeroPointYmm': +82.162, 'angle': pi},
    'BF' : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 2255, 'zeroPointXmm':  -1.0, 'zeroPointYmm': +82.162, 'angle': 3*pi/2},
    'CF' : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 2255, 'zeroPointXmm':  -1.0, 'zeroPointYmm': +82.162, 'angle': pi/2},
    'DF' : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 2255, 'zeroPointXmm':  -1.0, 'zeroPointYmm': +82.162, 'angle': 0}
}







def sunSkyCoordinates(julianDate):

    """
    PURPOSE: Compute the equatorial sky coordinates of the Sun at the given time.

    INPUT: julianDate: Julian date (floating point). The distinction between JD and BJD is neglected.
                       E.g. 2455000.5

    OUTPUT: (RA, Dec): equatorial sky coordinates of the Sun at the given time [rad]

    REMARK: Source: "Computing the solar vector", Blanco-Muriel et al., (2001), Solar Energy, Vol 70, pp 431-441
    """

    # Compute the (fractional) number of days since 1 Jan 2000 at 12:00 noon.

    elapsedJulianDays = julianDate - 2451545.0

    # In the following we assume that the solar ecliptic latitude is always exactly 0.0.

    Omega = 2.1429 - 0.0010394594 * elapsedJulianDays
    meanLongitude = 4.8950630 + 0.017202791698 * elapsedJulianDays
    meanAnomaly = 6.2400600 + 0.0172019699 * elapsedJulianDays

    eclipticLongitude = meanLongitude + 0.03341607 * sin(meanAnomaly) + 0.00034894 * sin(2*meanAnomaly) - 0.0001134 - 0.0000203 * sin(Omega)
    eclipticObliquity = 0.4090928 - 6.2140e-9 * elapsedJulianDays + 0.00003963 * cos(Omega)

    # Compute the RA, DEC of the Sun. Ensure that the RA is positive.

    rightAscensionSun = np.arctan2(cos(eclipticObliquity) * sin(eclipticLongitude), cos(eclipticLongitude))
    declinationSun = np.arcsin(sin(eclipticObliquity) * sin(eclipticLongitude))

    if rightAscensionSun < 0.0:
        rightAscensionSun += 2.0 * np.pi

    # That's it

    return rightAscensionSun, declinationSun











def ecliptic2equatorial(lam, beta):

    """
    Convert ecliptic coordinates (lambda, beta) into equatorial coordinates (alpha, delta)

    INPUT: lam:   ecliptic longitude        [rad]
           beta:  ecliptic latitude         [rad]
 
    OUTPUT: alpha: equtorial right ascension [rad]
            delta: equatorial declination    [rad]
    """

    obliquity = 0.409087723                  # Obliquity of the ecliptic = 23.439 deg  [rad]
 
    sindelta = sin(beta) * cos(obliquity) + cos(beta) * sin(obliquity) * sin(lam)
    delta = arcsin(sindelta)
    cosdelta = cos(delta)

    if (cosdelta == 0.0):
        print("ecliptic2equatorial: pointing to equatorial pole.")
        return None,None
    
    sinalpha = (-sin(beta) * sin(obliquity) + cos(beta) * cos(obliquity) * sin(lam)) / cosdelta
    cosalpha = cos(lam) * cos(beta) / cosdelta

    alpha = arctan2(sinalpha, cosalpha)

    if (alpha < 0.0): alpha += 2.0 * 3.141592653589793

    return alpha, delta















def equatorial2ecliptic(alpha, delta):

    """
    Convert equatorial coordinates (alpha, delta) into ecliptic coordinates (lambda, beta)

    INPUT: alpha: equatorial right ascension [rad]
           delta: equatorial declination [rad]

    OUTPUT: lam:  ecliptic longitude [rad]
            beta: ecliptic latitude [rad]
    """

    obliquity = 0.409087723                   # Obliquity of the ecliptic = 23.439 deg  [rad]
 
    sinbeta = sin(delta) * cos(obliquity) - cos(delta) * sin(obliquity) * sin(alpha)
    beta = arcsin(sinbeta)
    cosbeta = cos(beta)

    if (cosbeta == 0.0):
        print("equatorial2ecliptic: pointing to ecliptic pole.")
        return None, None

    sinlambda = (sin(delta) * sin(obliquity) + cos(delta) * cos(obliquity) * sin(alpha)) / cosbeta
    coslambda = cos(alpha) * cos(delta) / cosbeta

    lam = arctan2(sinlambda, coslambda)

    if (lam < 0.0): lam += 2.0 * 3.141592653589793

    return lam, beta











def sunSkyCoordinatesAwayfromPlatformPointing(raPlatform, decPlatform):

    """
    Derive the location of the Sun which we assume to always be 180 degrees away from the platform pointing
    in the middle of the total time series.

    INPUT: raPlatform:  right ascension of the pointing of the Platform [rad]
           decPlatform: declination of the pointing of the Platform [rad]

    OUTPUT: raSun:  right ascension of the sun [rad]
            decSun: declination of the sun [rad]
    """

    lambdaPlatform, betaPlatform = equatorial2ecliptic(raPlatform, decPlatform)
    lambdaSun = lambdaPlatform - np.pi
    if (lambdaSun < 0.0): lambdaSun += 2.0 * np.pi
    raSun, decSun = ecliptic2equatorial(lambdaSun, 0.0)

    return raSun, decSun











def skyToFocalPlaneCoordinates(raStar, decStar, raPlatform, decPlatform, tiltAngle, azimuthAngle, focalPlaneAngle, focalLength):

    """
    PURPOSE: Convert the equatorial sky coordinates (alpha, delta) of a star to undistorted normalized focal plane coordinates (xFP, yFP),
             assuming a spherical pinhole camera.

    INPUT: raStar:          right ascension of the star                               [rad]
           decStar:         declination of the star                                   [rad]
           raPlatform:      right ascension of the optical axis                       [rad]
           decPlatform:     declination of the optical axis                           [rad]
           tiltAngle:       tilt angle of the telescope w.r.t. platform z-axis        [rad]                   
           azimuthAngle:    azimuth angle of the telescope on the platform            [rad]
           focalPlaneAngle: angle between the Y_TL axis and the Y_FP axis: gamma_FP   [rad]
           focalLength:     focal length of the camera. Unit: see remarks.

    OUTPUT: xFP, yFP: normalized cartesian coordinates of the project star in the focal plane in the FP reference frame. 
                      Unit: see remarks:
    
    REMARK: - The unit of the cartesian focal plane coordinates is the same as the one of focalLength. 
              focalLength can be expressed in e.g. [mm] or in [pixels]. If focalLength == 1.0, then the corresponding
              focal plane coordinates are called "normalized coordinates."
            - Reference documents: PLATO-KUL-PL-TN-0001 and PLATO-DLR-PL-TN-016

    """

    # Get the sky position of the Sun (ra, dec) [rad]

    raSun, decSun = sunSkyCoordinatesAwayfromPlatformPointing(raPlatform, decPlatform)

    # Compute the equatorial cartesian coordinates of the unit vector along the z-axis (= roll = pointing axis) of the platform.
    # The x-axis of the platform points to the highest point fof the sunshield, which is pointing to the (average) sky position
    # of the Sun.

    zSC = np.array([cos(decPlatform)*cos(raPlatform), cos(decPlatform)*sin(raPlatform), sin(decPlatform)])
    deltax = np.arctan(- cos(raPlatform-raSun) / tan(decPlatform))
    xSC = np.array([cos(deltax)*cos(raSun), cos(deltax)*sin(raSun), sin(deltax)])
    ySC = np.cross(zSC, xSC)

    # Compute the rotation matrix to convert cartesian coordinates in the equatorial reference frame to 
    # cartesian coordinates in the spacecraft framework.

    rotEQ2SC   = np.array([[xSC[0], xSC[1], xSC[2]], \
                           [ySC[0], ySC[1], ySC[2]], \
                           [zSC[0], zSC[1], zSC[2]]])

    
    # Compute the rotation matrix to convert cartesian coordinates in the spacecraft reference frame to 
    # cartesian coordinates in the telescope reference frame

    rotAzimuth = np.array([[cos(azimuthAngle), -sin(azimuthAngle), 0],   \
                           [sin(azimuthAngle),  cos(azimuthAngle), 0],   \
                           [        0        ,          0,         1]])

    rotTilt = np.array([[cos(tiltAngle), 0, -sin(tiltAngle)], \
                        [     0        , 1,        0       ], \
                        [sin(tiltAngle), 0,  cos(tiltAngle)]])

    rotSC2TL = np.dot(rotTilt, rotAzimuth)

    # Compute the rotation matrix to convert cartesian coordinates in the telescope reference frame to
    # cartesian coordinates in the focal plane reference frame

    rotTL2FP = np.array([[ cos(focalPlaneAngle), sin(focalPlaneAngle), 0],  \
                         [-sin(focalPlaneAngle), cos(focalPlaneAngle), 0],  \
                         [          0          ,           0         , 1]])

    # Combine all the rotation matrices

    rotEQ2FP = np.dot(rotTL2FP, np.dot(rotSC2TL, rotEQ2SC))

    # Compute the cartesian coordinates of the star in the equatorial reference frame

    starEQ = np.array([cos(decStar)*cos(raStar), cos(decStar)*sin(raStar), sin(decStar)])

    # Transform these coordinates to the corresponding ones in the focal plane reference frame:

    starFP = np.dot(rotEQ2FP, starEQ)

    # Convert the units to the one of focalLength (usually [mm]), and normalize the coordinates 
    # to take into account the pinhole camera projection.

    xFPmm = - focalLength * starFP[0]/starFP[2]
    yFPmm = - focalLength * starFP[1]/starFP[2]

    # That's it

    return xFPmm, yFPmm
    











def focalPlaneToSkyCoordinates(xFP, yFP, raPlatform, decPlatform, tiltAngle, azimuthAngle, focalPlaneAngle, focalLength):

    """
    PURPOSE: Convert the undistorted normalized focal plane coordinates (xFP, yFP) of a star to equatorial sky coordinates (alpha, delta),
             assuming a spherical pinhole camera.

    INPUT: xFP:             undistorted normalized cartesian x-coordinate in the focal plane reference frame [same unit as focalLength]
           yFP:             undistorted normalized cartesian y-coordinate in the focal plane reference frame [same unit as focalLength]
           raPlatform:      right ascension of the platform pointing axis             [rad]
           decPlatform:     declination of the platform pointing axis                 [rad]
           tiltAngle:       tilt angle of the telescope w.r.t. platform z-axis        [rad]                   
           azimuthAngle:    azimuth angle of the telescope on the platform            [rad]
           focalPlaneAngle: angle between the Y_TL axis and the Y_FP axis: gamma_FP   [rad]
           focalLength:     focal length of the camera.                               [mm]

    OUTPUT: raStar, decStar: Equatorial sky coordinates, right ascension and declination, of the star [rad]

    REMARK: The transformation assumes that the pinhole reverses the image.
    """

    # Get the sky position of the Sun (ra, dec) [rad]

    raSun, decSun = sunSkyCoordinatesAwayfromPlatformPointing(radPlatform, decPlatform)

    # Undo the reverse-image projection effect of the pinhole

    vecFP = array([-xFP/focalLength, -yFP/focalLength, 1.0])

    # Compute the rotation matrix to convert cartesian coordinates in the focal plane reference frame to
    # cartesian coordinates in the telescope reference frame

    rotFP2TL = np.array([[cos(focalPlaneAngle), -sin(focalPlaneAngle), 0],  \
                         [sin(focalPlaneAngle),  cos(focalPlaneAngle), 0],  \
                         [          0          ,           0         , 1]])

    # Compute the rotation matrix to convert cartesian coordinates in the telescope reference frame to 
    # cartesian coordinates in the spacecraft reference frame

    rotAzimuth = np.array([[ cos(azimuthAngle), sin(azimuthAngle), 0],   \
                           [-sin(azimuthAngle), cos(azimuthAngle), 0],   \
                           [        0        ,          0,         1]])

    rotTilt = np.array([[ cos(tiltAngle), 0, sin(tiltAngle)], \
                        [     0        , 1,        0       ], \
                        [-sin(tiltAngle), 0, cos(tiltAngle)]])

    rotTL2SC = np.dot(rotAzimuth, rotTilt)


    # Compute the equatorial cartesian coordinates of the unit vector along the z-axis (= roll = pointing axis) of the platform.
    # The x-axis of the platform points to the highest point fof the sunshield, which is pointing to the (average) sky position
    # of the Sun.

    zSC = np.array([cos(decPlatform)*cos(raPlatform), cos(decPlatform)*sin(raPlatform), sin(decPlatform)])
    deltax = np.arctan(- cos(raPlatform-raSun) / tan(decPlatform))
    xSC = np.array([cos(deltax)*cos(raSun), cos(deltax)*sin(raSun), sin(deltax)])
    ySC = np.cross(zSC, xSC)

    # Compute the rotation matrix to convert cartesian coordinates in the equatorial reference frame to 
    # cartesian coordinates in the spacecraft reference frame

    rotSC2EQ   = np.array([[xSC[0], ySC[0], zSC[0]], \
                           [xSC[1], ySC[1], zSC[1]], \
                           [xSC[2], ySC[2], zSC[2]]])
    
    # Combine all the rotation matrices

    rotFP2EQ = np.dot(rotSC2EQ, np.dot(rotTL2SC, rotFP2TL))

    # Transform the unnormalized focal plane coordinates to the corresponding ones in the equatorial reference frame

    vecEQ = np.dot(rotFP2EQ, vecFP)

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
















def undistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm, distortionCoefficients):
    
    """
    PURPOSE:      Convert from undistorted to distorted normalized focal plane coordinates
    
    INPUTS:       xFPmm  undistorted normalized focal plane x-coordinate [mm]
                  yFPmm  undistorted normalized focal plane y-coordinate [mm]
                  distortionCoefficients  List of polynomial coefficients
    
    OUTPUTS:      (xFPdist, yFPdist) distorted x and y coordinates [mm]

    Note: Example of distortion coefficients: [-0.0036696919678, 1.0008542317, -4.12553764817e-05, 5.7201219949e-06]
    """

    P = Polynomial(distortionCoefficients)

    angle = arctan2(yFPmm, xFPmm)    # [radians]
    
    rFP = sqrt(xFPmm * xFPmm + yFPmm * yFPmm)
    rFPdist = P(rFP)
    
    xFPdist = cos(angle) * rFPdist
    yFPdist = sin(angle) * rFPdist
    
    return xFPdist, yFPdist










def distortedToUndistortedFocalPlaneCoordinates(xFPdist, yFPdist, inverseDistortionCoefficients):
    
    """
    PURPOSE:     Convert from distorted to undistorted normalized focal plane coordinates
    
    INPUTS:      xFPdist  Distorted normalized focal plane x-coordinate [mm]
                 yFPdist  DIstorted normalized focal plane y-coordinate [mm]
                 inverseDistortionCoefficients  List of polynomial coefficients
    
    OUTPUTS:     (xFPmm, yFPmm) distorted x and y coordinates [mm]

    Note: Example of inverse distortion coefficients: [-0.00458067036444, 1.00110311283, -5.61136295937e-05, -4.311925329e-06]
    """
    
    IP = Polynomial(inverseDistortionCoefficients)
    
    angle = arctan2(yFPdist, xFPdist)  # [radians]
    
    rFP   = sqrt(xFPdist * xFPdist + yFPdist * yFPdist)
    rFPmm = IP(rFP)
    
    xFPmm = cos(angle) * rFPmm
    yFPmm = sin(angle) * rFPmm
    
    return xFPmm, yFPmm











def pixelToFocalPlaneCoordinates(xCCDpixel, yCCDpixel, pixelSize, ccdZeroPointX, ccdZeroPointY, CCDangle):

    """
    PUROSE: Given the (real-valued) pixel coordinates of the star on the CCD, compute the (x,y)
            coordinates in the FP reference system.

    INPUT: xCCDpixel     : x-coordinate (column-number !!) of the star on the CCD  [pixel]
           yCCDpixel     : y-coordinate (row-number !!)  of the star on the CCD  [pixel]
           pixelSize     : size of 1 pixel in micron (not [mm]!)
           ccdZeroPointX : x-coordinate of the CCD (0,0) point in the FP' reference system [mm]
           ccdZeroPointY : y-coordinate of the CCD (0,0) point in the FP' reference system [mm]
           CCDangle      : CCD orientation angle in the FP' reference frame  [rad]

    OUTPUT: xFP: column pixel coordinate of the star (real-valued) [mm]
            yFP: row pixel coordinate of the star (real-valued) [mm]
    """

    # Convert the pixel coordinates into [mm] coordinates

    xCCDmm = xCCDpixel * pixelSize / 1000.0
    yCCDmm = yCCDpixel * pixelSize / 1000.0

    # Convert the CCD coordinates into FP coordinates [mm]

    xFP = (xCCDmm - ccdZeroPointX) * cos(CCDangle) - (yCCDmm - ccdZeroPointY) * sin(CCDangle)
    yFP = (xCCDmm - ccdZeroPointX) * sin(CCDangle) + (yCCDmm - ccdZeroPointY) * cos(CCDangle)

    # That's it

    return xFP, yFP











def focalPlaneToPixelCoordinates(xFP, yFP, pixelSize, ccdZeroPointX, ccdZeroPointY, CCDangle):

    """
    PUROSE: Compute the (real-valued) pixel coordinates of the star on the CCD, given the (x,y)
            coordinates in the FP reference system. The cartesian focal plane coordinates are
            supposed to be normalized (i.e. with the pinhole projection taken into account).

    INPUT: xFP           : normalized x-coordinate of the star in the FP reference system  [mm]
           yFP           : normalized y-coordinate of the star in the FP reference system  [mm]
           pixelSize     : size of 1 pixel [micron]  (not [mm]!)
           ccdZeroPointX : x-coordinate of the CCD (0,0) point in the FP reference system [mm]
           ccdZeroPointY : y-coordinate of the CCD (0,0) point in the FP reference system [mm]
           CCDangle      : CCD orientation angle in the FP reference frame  [rad]

    OUTPUT: xCCDpixel: column pixel coordinate of the star (real-valued) 
            yCCDpixel: row pixel coordinate of the star (real-valued)
    """

    # Convert the FP coordinates into CCD coordinates [mm]

    xCCDmm = ccdZeroPointX + xFP * cos(CCDangle) + yFP * sin(CCDangle)
    yCCDmm = ccdZeroPointY - xFP * sin(CCDangle) + yFP * cos(CCDangle)

    # Convert the [mm] coordinates into pixel coordinates

    xCCDpixel = xCCDmm / pixelSize * 1000.0
    yCCDpixel = yCCDmm / pixelSize * 1000.0

    # That's it

    return xCCDpixel, yCCDpixel













def gnomonicRadialDistanceFromOpticalAxis(xFP, yFP, focalLength):

    """
    Calculate the gnomonic radial distance with respect to the optical axis in the focal plane

    INPUT: xFP  Focal plane x-coordinate [mm]
           yFP  Focal plane y-coordinate [mm]
 
    OUTPUT: the angular distance of the star w.r.t. the optical axis [rad]

    """

    tanx = xFP / focalLength
    tany = yFP / focalLength

    angularDistance = arccos(1.0/sqrt(1.0 + tanx*tanx + tany*tany));

    # Take care that the angle is between [0, 2*PI]

    if angularDistance < 0.0:
        angularDistance += 2.0 * np.pi
    elif angularDistance > 2.0 * np.pi:
        angularDistance -= 2.0 * np.pi

    # That's it!

    return angularDistance;












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












def getCCDandPixelCoordinates(raStar, decStar, raPlatform, decPlatform, tiltAngle, azimuthAngle,  \
                              focalPlaneAngle, focalLength, pixelSize, includeFieldDistortion, distortionCoefficients, \
                              normal):

    """
    PURPOSE: Given the equatorial coordinates of a star, find out on which CCD it falls ('A', 'B', ...)
             and compute the pixel coordinates of the star on this CCD. If the star doesn't fall on any of the CCDs
             then (None, None, None) is given as output.

    INPUT: raStar:          right ascension of the star                               [rad]
           decStar:         declination of the star                                   [rad]      
           raPlatform:      right ascension of the platform pointing axis             [rad]
           decPlatform:     declination of the platform pointing axis                 [rad]
           tiltAngle:       tilt angle of the telescope w.r.t. platform z-axis        [rad]                   
           azimuthAngle:    azimuth angle of the telescope on the platform            [rad]
           focalPlaneAngle: angle between the Y_TL axis and the Y_FP axis: gamma_FP   [rad]
           focalLength:     focal length of the camera.                               [mm]
           pixelSize:       pixel size                                                [micron]
           includeFieldDistortion: True to include field distortion in coordinate transformations, false otherwise
           distortionCoefficients: Coefficients of the polynomial describing the distortion 
           normal:                 True for the normal camera configuration, False for the fast cameras

    OUTPUT: ccdCode: for normal camera: either 'A', 'B', 'C', or 'D'
                     for fast camer: either 'AF', 'BF', 'CF', 'DF'
                     if on no CCD: None    
            xCCDpix: x-coordinate (column number) of the star on the CCD  [pix]
                     if on no CCD: None
            yCCDpix: y-coordinate (row number) of the star on the CCD  [pix]
                     if on no CCD: None 

    Note: Example of distortion coefficients: [-0.0036696919678, 1.0008542317, -4.12553764817e-05, 5.7201219949e-06]            
    """

    # Select the proper CCD codes depending on whether we're dealing with the nominal or the fast cams

    if normal == True:
        ccdCodes = ['A', 'B', 'C', 'D']
    else:
        ccdCodes = ['AF', 'BF', 'CF', 'DF']


    # Compute the (x,y) coordinates in the FP reference system [mm]

    xFPmm, yFPmm = skyToFocalPlaneCoordinates(raStar, decStar, raPlatform, decPlatform, tiltAngle, azimuthAngle, focalPlaneAngle, focalLength)

    if (includeFieldDistortion == True) or (includeFieldDistortion == "yes"):
        xFPmm, yFPmm = undistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm, distortionCoefficients)

    # Concert to CCD pixel coordinates
    # Note that the pixel size is in [micron] not in [mm]

    xCCDpixel = xFPmm / pixelSize * 1000.0
    yCCDpixel = yFPmm / pixelSize * 1000.0

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












def platformToTelescopePointingCoordinates(raPlatform, decPlatform, raSun, decSun, azimuthAngle, tiltAngle):

    """
    PURPOSE: Given the platform pointing coordinates (i.e. the sky coordinates of the jitter axis)
             and the orientation of the telescope on the platform, compute the sky coordinates of
             the optical axis of the telescope.
             See also: PLATO-KUL-PL-TN-001

    INPUT: raPlatform:   Right Ascension of the Platform pointing axis           [rad]
           declatform:   Declination of the platfor pointing axis                [rad]
           raSun:        Right Ascension of the Sun (shield)                     [rad]
           decSun:       Declination of the Sun (shield)                         [rad]
           azimuthAngle: Azimuth of the telescope on the platform                [rad]
           tiltAngle:    Tilt angle between platform and telescope pointing axes [rad]

    OUTPUT: raTelescope: right ascension of the optical axis of the telescope [rad]
            decTelescope: declination of the optical axis of the telescope [rad]

    """


    # Compute the rotation matrix to convert cartesian coordinates in the telescope reference frame to 
    # cartesian coordinates in the spacecraft reference frame

    rotAzimuth = np.array([[ cos(azimuthAngle), sin(azimuthAngle), 0],   \
                           [-sin(azimuthAngle), cos(azimuthAngle), 0],   \
                           [        0        ,          0,         1]])

    rotTilt = np.array([[ cos(tiltAngle), 0, sin(tiltAngle)], \
                        [     0        , 1,        0       ], \
                        [-sin(tiltAngle), 0, cos(tiltAngle)]])

    rotTL2SC = np.dot(rotAzimuth, rotTilt)


    # Compute the equatorial cartesian coordinates of the unit vector along the z-axis (= roll = pointing axis) of the platform.
    # The x-axis of the platform points to the highest point fof the sunshield, which is pointing to the (average) sky position
    # of the Sun.

    zSC = np.array([cos(decPlatform)*cos(raPlatform), cos(decPlatform)*sin(raPlatform), sin(decPlatform)])
    deltax = np.arctan(- cos(raPlatform-raSun) / tan(decPlatform))
    xSC = np.array([cos(deltax)*cos(raSun), cos(deltax)*sin(raSun), sin(deltax)])
    ySC = np.cross(zSC, xSC)

    # Compute the rotation matrix to convert cartesian coordinates in the equatorial reference frame to 
    # cartesian coordinates in the spacecraft reference frame

    rotSC2EQ   = np.array([[xSC[0], ySC[0], zSC[0]], \
                           [xSC[1], ySC[1], zSC[1]], \
                           [xSC[2], ySC[2], zSC[2]]])
    
    # Combine all the rotation matrices

    rotTL2EQ = np.dot(rotSC2EQ, rotTL2SC)

    # In the telescope reference frame, the optical axis is simply the z-axis = (0,0,1)

    vecTL = np.array([0.0, 0.0, 1.0])

    # Get the equatorial coordinates of this optical axis vector.

    vecEQ = np.dot(rotTL2EQ, vecTL)

    # Convert the cartesian equatorial coordinates to equatorial sky coordinates

    norm = sqrt(vecEQ[0]*vecEQ[0] + vecEQ[1]*vecEQ[1] + vecEQ[2]*vecEQ[2]) 
    decOpticalAxis = pi/2.0 - arccos(vecEQ[2]/norm)
    raOpticalAxis = arctan2(vecEQ[1], vecEQ[0])

    # Ensure that the right ascension is positive

    if isinstance(raOpticalAxis, np.ndarray):
        raOpticalAxis[raOpticalAxis < 0.0] += 2.*pi
    else:
        if (raOpticalAxis < 0.0):
            raOpticalAxis += 2.*np.pi

    # That's it

    return (raOpticalAxis, decOpticalAxis)












def calculateSubfieldAroundCoordinates(subfieldSizeX, subfieldSizeY, raStar, decStar, raPlatform, decPlatform, \
                                       tiltTelescope, azimuthTelescope, focalPlaneAngle, focalLength, pixelSize,              \
                                       includeFieldDistortion, distortionCoefficients, normal):

    """
    PURPOSE: Calculates the location of the subfield such that the star with coordinates (raStar, decStar)
             is centered in the subfield. The function also checks if there is enough space around the pixel
             to accommodate a subfield of the specified size. 

    NOTE:    This function is used by setSubfieldAroundCoordinates() and usually does not need to be called by the user. 

    INPUT: subfieldSizeX:          full width (# of columns) of the subfield                [pix]
           subfieldSizeY:          full height (#of rows) of the subfield                   [pix]
           raStar:                 right ascension of the star                              [rad]
           decStar:                declination of the star                                  [rad]
           raPlatform:             right ascension of the platform pointing axis            [rad]
           decPlatform:            declination of the platform pointing axis                [rad]
           tiltTelescope           tilt angle of the telescope w.r.t. platform z-axis       [rad]                   
           azimuthTelescope:       azimuth angle of the telescope on the platform           [rad]
           focalPlaneAngle:        angle between the Y_TL axis and the Y_FP axis: gamma_FP  [rad]
           focalLength:            focal length of the camera.                              [mm]
           pixelSize:              pixel size                                               [micron]
           includeFieldDistortion: True to include field distortion in coordinate transformations, false otherwise
           distortionCoefficients: Coefficients of the polynomial describing the distortion
           normal:                 True for the normal camera configuration, False for the fast cameras

    OUTPUT: ccdCode: "A", "B", "C" or "D" if nominal=True, "AF", "BF", "CF" or "DF" otherwise
            xCCDpix: x-coordinate of the star in pixels (i.e. column number)
            yCCDpix: y-coordinate of the star in pixels (i.e. row number)

    REMARKS: - If the coordinates do not fall on any CCD, an error message is shown
             - If the star is too close to the edge for the given subfield size, and error message is shown,
               followed by an exit(1)
             - Example of distortion coefficients: [-0.0036696919678, 1.0008542317, -4.12553764817e-05, 5.7201219949e-06]
    """

    # Find out on which CCD the star falls, and the corresponding pixel coordinates

    ccdCode, xCCDpix, yCCDpix = getCCDandPixelCoordinates(raStar, decStar, raPlatform, decPlatform, tiltTelescope, \
                                                          azimuthTelescope, focalPlaneAngle, focalLength, pixelSize, \
                                                          includeFieldDistortion, distortionCoefficients, normal)

    # If the CCD code is None, the star does not fall on any ccd -> error

    if ccdCode == None:
        print("Error: coordinates do not fall on any CCD")
        print("raStar, decStar = {0}, {1}".format(raStar, decStar))
        return None, None, None


    # If the star does fall on a CCD, check if it's not too close to the edge for the subfield to
    # be completely on the CCD.

    xCCDpix = round(xCCDpix)                # integer values
    yCCDpix = round(yCCDpix)
    firstRow = CCD[ccdCode]["firstRow"]     # different from nominal than for fast cams
    Ncols = CCD[ccdCode]["Ncols"]
    Nrows = CCD[ccdCode]["Nrows"]

    if     (xCCDpix - subfieldSizeX/2 < 0)        or (xCCDpix + subfieldSizeX/2 - 1 > Ncols-1)   \
        or (yCCDpix - subfieldSizeY/2 < firstRow) or (yCCDpix + subfieldSizeY/2 - 1 > Nrows-1): 

        print("Error: pixel coordinates (row, col) = ({0},{1}) too close to the edge to accommodate subfield with size {2}x{3}" \
           .format(yCCDpix,xCCDpix,subfieldSizeX, subfieldSizeY))
        return None, None, None

    # That's it!

    return ccdCode, xCCDpix, yCCDpix












def skyToPixelCoordinates(sim, raStar, decStar, normal):
    """
    PURPOSE: Convert sky coordinates to pixel coordinates
    
    NOTE:   It is assumed that the configuration parameters in the sim object contains
            a correct (ra, dec) of the platform, a correct (azimuth, tilt) of the telescope,
            a valid value for the focal length, the plate scale, the pixel size, and that
            the switch to include distortion or not is set correctly.
            
    INPUT:  sim:            instance of simulation class, see simulation.py
            raStar:                 right ascension of the star                    [rad]
            decStar:                declination                                    [rad]
            subfieldSizeX:          width (i.e. number of columns) of the subiield [pix]
            subfieldSizeY:          height (i.e. number of rows) of the sub-field  [pix]
            normal:                 True for the normal camera configuration, False for the fast cameras
    
    OUTPUT: ccdCode
            xCCDpixel: column pixel coordinate of the star (real-valued) 
            yCCDpixel: row pixel coordinate of the star (real-valued)
    """
    
    if (sim["Camera/IncludeFieldDistortion"] == "yes")  or (sim["Camera/IncludeFieldDistortion"] == "1"):
        includeFieldDistortion = True
        distortionCoefficients = sim["Camera/FieldDistortion/Coefficients"]
    else:
        includeFieldDistortion = False

    pixelSize        = float(sim["CCD/PixelSize"])
    focalLength      = float(sim["Camera/FocalLength"]) * 1000.0                   # [m] -> [mm]
    raPlatform       = np.deg2rad(float(sim["ObservingParameters/RApointing"]))
    decPlatform      = np.deg2rad(float(sim["ObservingParameters/DecPointing"]))
    focalPlaneAngle  = float(sim["Camera/FocalPlaneOrientation"])
    azimuthTelescope = np.deg2rad(float(sim["Telescope/AzimuthAngle"]))
    tiltTelescope    = np.deg2rad(float(sim["Telescope/TiltAngle"]))
    
    # Get the sky position of the Sun (ra, dec) in the middle of the time series [rad] 

    raSun, decSun = sunSkyCoordinatesAwayfromPlatformPointing(radPlatform, decPlatform)

    # Get the pixel coordinates on the CCD

    ccdCode, xCCDpixel, yCCDpixel = getCCDandPixelCoordinates(raStar, decStar, raSun, decSun, raPlatform, decPlatform, tiltTelescope, azimuthTelescope,  \
                                                              focalPlaneAngle, focalLength, pixelSize, includeFieldDistortion, distortionCoefficients, normal)

    return ccdCode, xCCDpixel, yCCDpixel
















def pixelToSkyCoordinates(sim, ccdCode, xCCDpixel, yCCDpixel):

    """
    PURPOSE: Convert pixel coordinates to equatorial sky coordinates (ra,dec)
    
    NOTE:   It is assumed that the configuration parameters in the sim object contains
            a correct (ra, dec) of the platform, a correct (azimuth, tilt) of the telescope,
            a valid value for the focal length, the pixel size, and that
            the switch to include distortion or not is set correctly.
            
    INPUT:  sim:        simulation for which the configuration file is adapted
            ccdCode:    for nominal camera: either 'A', 'B', 'C', 'D'
                        for fast camera: either 'AF', 'BF', 'CF', 'DF'
            xCCDpixel:  x-coordinate (column-number) of the star on the CCD  [pix]
            yCCDpixel:  y-coordinate (row-number) of the star on the CCD     [pix]
    
    OUTPUT: raStar, decStar: Equatorial coordinates (right ascension and declination) of the star [rad]
    """
    
    if (sim["Camera/IncludeFieldDistortion"] == "yes")  or (sim["Camera/IncludeFieldDistortion"] == "1"):
        includeFieldDistortion = True
        inverseDistortionCoefficients = sim["Camera/FieldDistortion/InverseCoefficients"]
    else:
        includeFieldDistortion = False

    pixelSize        = float(sim["CCD/PixelSize"])
    focalLength      = float(sim["Camera/FocalLength"]) * 1000.0                     # [m] -> [mm]
    raPlatform       = np.deg2rad(float(sim["ObservingParameters/RApointing"]))
    decPlatform      = np.deg2rad(float(sim["ObservingParameters/DecPointing"]))
    focalPlaneAngle  = float(sim["Camera/FocalPlaneOrientation"])
    azimuthTelescope = np.deg2rad(float(sim["Telescope/AzimuthAngle"]))
    tiltTelescope    = np.deg2rad(float(sim["Telescope/TiltAngle"]))

    ccdZeroPointX = CCD[ccdCode]['zeroPointXmm']
    ccdZeroPointY = CCD[ccdCode]['zeroPointYmm']
    ccdAngle      = CCD[ccdCode]['angle']
    

    # Get the sky position of the Sun (ra, dec) in the middle of the time series [rad] 

    raSun, decSun = sunSkyCoordinatesAwayfromPlatformPointing(radPlatform, decPlatform)

    # Get the focal plane coordinates

    xFPmm, yFPmm = pixelToFocalPlaneCoordinates(xCCDpixel, yCCDpixel, pixelSize, ccdZeroPointX, ccdZeroPointY, ccdAngle)
    
    # If required, undistort them
    
    if includeFieldDistortion:
        xFPmm, yFPmm = distortedToUndistortedFocalPlaneCoordinates(xFPmm, yFPmm, inverseDistortionCoefficients)
    
    # Get the corresponding sky coordinats

    ra, dec = focalPlaneToSkyCoordinates(xFPmm, yFPmm, raSun, decSun, raPlatform, decPlatform, tiltTelescope, azimuthTelescope, \
                                         focalPlaneAngle, focalLength)

    return ra, dec
