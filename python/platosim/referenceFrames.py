#!/usr/bin/env python3

import h5py
import math
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy.polynomial import Polynomial
from numpy.linalg import norm


# CCD configuration
#
# ccdCode:      1, 2, 3, 4: nominal cameras; 1F, 2F, 3F, 4F: fast cameras
# Ncols:        Number of exposed columns (column number varies along x-coordinate)
# Nrows:        Number of exposed rows (row number varies along y-coordinate)
# firstRow:     First row that is exposed. For the nominal cams this is simply row 0.
#               For the fast cams, the exposed rows are rows 2255 until 4510, after the exposure
#               are then frame-transfered to rows 0 until 2254.
# zeroPointXmm: x-coordinate of the (0,0) pixel of the CCD, in the FP' reference frame [mm]
# zeroPointYmm: y-coordinate of the (0,0) pixel of the CCD, in the FP' reference frame [mm]
# angle:        gamma_{ccd} [rad]: orientation angle of the CCD in the FP' reference frame

CCD = \
{
    '1'  : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 0,    'zeroPointXmm':  -1.3, 'zeroPointYmm': +82.48, 'angle': np.pi},
    '2'  : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 0,    'zeroPointXmm':  -1.3, 'zeroPointYmm': +82.48, 'angle': 3*np.pi/2},
    '3'  : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 0,    'zeroPointXmm':  -1.3, 'zeroPointYmm': +82.48, 'angle': 0},
    '4'  : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 0,    'zeroPointXmm':  -1.3, 'zeroPointYmm': +82.48, 'angle': np.pi/2},
    '1F' : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 2255, 'zeroPointXmm':  -1.3, 'zeroPointYmm': +82.48, 'angle': np.pi},
    '2F' : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 2255, 'zeroPointXmm':  -1.3, 'zeroPointYmm': +82.48, 'angle': 3*np.pi/2},
    '3F' : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 2255, 'zeroPointXmm':  -1.3, 'zeroPointYmm': +82.48, 'angle': 0},
    '4F' : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 2255, 'zeroPointXmm':  -1.3, 'zeroPointYmm': +82.48, 'angle': np.pi/2}
}






def sunSkyCoordinates(julianDate):

    """Compute the equatorial sky coordinates of the Sun at the given time.

    Source
    ------
    "Computing the solar vector", 
    Blanco-Muriel et al. (2001), Solar Energy, Vol 70, pp 431-441

    Parameters
    ----------
    julianDate : float
        Julian date (e.g. 2455000.5). The distinction between JD and BJD is neglected.
                      
    Return
    ------
    raSun : float
        Right ascension of the Sun at the given time [rad]
    decSun : float
        Declination of the Sun at the given time [rad]
    """

    # Compute the (fractional) number of days since 1 Jan 2000 at 12:00 noon.

    elapsedJulianDays = julianDate - 2451545.0

    # In the following we assume that the solar ecliptic latitude is always exactly 0.0.

    Omega         = 2.1429    - 0.0010394594   * elapsedJulianDays
    meanLongitude = 4.8950630 + 0.017202791698 * elapsedJulianDays
    meanAnomaly   = 6.2400600 + 0.0172019699   * elapsedJulianDays

    eclipticLongitude = (meanLongitude + 0.03341607 * np.sin(meanAnomaly) + 0.00034894 *
                         np.sin(2*meanAnomaly) - 0.0001134 - 0.0000203 * np.sin(Omega))
    eclipticObliquity = 0.4090928 - 6.2140e-9 * elapsedJulianDays + 0.00003963 * np.cos(Omega)

    # Compute the RA, DEC of the Sun. Ensure that the RA is positive.

    raSun  = np.arctan2(np.cos(eclipticObliquity) * np.sin(eclipticLongitude), np.cos(eclipticLongitude))
    decSun = np.arcsin(np.sin(eclipticObliquity)  * np.sin(eclipticLongitude))

    if raSun < 0.0: raSun += 2.0 * np.pi

    # That's it

    return raSun, decSun





def ecliptic2equatorial(lam, beta):

    """Convert ecliptic coordinates (lambda, beta) into equatorial coordinates (alpha, delta)

    Parameters
    ----------
    lam : float
        Ecliptic longitude [rad]
    beta : float  
        Eecliptic latitude [rad]

    Return
    ------
    alpha : float
        Equtorial right ascension [rad]
    delta : float
        Equatorial declination [rad]
    """

    # Obliquity of the ecliptic = 23.439 deg in [rad]
    
    obliquity = 0.409087723

    sindelta = np.sin(beta) * np.cos(obliquity) + np.cos(beta) * np.sin(obliquity) * np.sin(lam)
    delta    = np.arcsin(sindelta)
    cosdelta = np.cos(delta)

    if (cosdelta == 0.0):
        print("ecliptic2equatorial: pointing to equatorial pole.")
        return None,None

    sinalpha = (-np.sin(beta) * np.sin(obliquity) + np.cos(beta) * np.cos(obliquity) * np.sin(lam)) / cosdelta
    cosalpha = np.cos(lam) * np.cos(beta) / cosdelta

    alpha = np.arctan2(sinalpha, cosalpha)

    if (alpha < 0.0): alpha += 2.0 * 3.141592653589793

    return alpha, delta





def equatorial2ecliptic(alpha, delta):

    """Convert equatorial coordinates (alpha, delta) into ecliptic coordinates (lambda, beta)

    Parameters
    ----------
    alpha : float
        Equtorial right ascension [rad]
    delta : float
        Equatorial declination [rad]

    Return
    ------
    lam : float
        Ecliptic longitude [rad]
    beta : float  
        Eecliptic latitude [rad]
    """

    # Obliquity of the ecliptic = 23.439 deg  [rad]
    
    obliquity = 0.409087723

    sinbeta = np.sin(delta) * np.cos(obliquity) - np.cos(delta) * np.sin(obliquity) * np.sin(alpha)
    beta    = np.arcsin(sinbeta)
    cosbeta = np.cos(beta)

    if (cosbeta == 0.0):
        print("equatorial2ecliptic: pointing to ecliptic pole.")
        return None, None

    sinlambda = (np.sin(delta) * np.sin(obliquity) + np.cos(delta) * np.cos(obliquity) * np.sin(alpha)) / cosbeta
    coslambda = np.cos(alpha) * np.cos(delta) / cosbeta

    lam = np.arctan2(sinlambda, coslambda)

    if (lam < 0.0): lam += 2.0 * 3.141592653589793

    return lam, beta





def sunSkyCoordinatesAwayfromPlatformPointing(raPlatform, decPlatform, solarPanelOrientation=0.0):

    """Location of Sun w.r.t. platform pointing.

    Derive the location of the Sun which we assume to always be 180 degrees
    away from the platform pointing in the middle of the total time series.

    Parameters
    ----------
    raPlatform : float
        Right ascension of the pointing of the Platform [rad]
    decPlatform : float 
        Declination of the pointing of the Platform [rad]
    solarPanelOrientation : float
        Orientation of the solar panel [rad]
        This corresponds to (0, pi/2, pi, 3pi/2) for quarters (Q1,Q2,Q3,Q4)

    Return
    ------
    raSun : float
        Right ascension of the Sun at the given time [rad]
    decSun : float
        Declination of the Sun at the given time [rad]
    """

    lambdaPlatform, betaPlatform = equatorial2ecliptic(raPlatform, decPlatform)
    lambdaSun = lambdaPlatform - np.pi + solarPanelOrientation

    if (lambdaSun < 0.0):
        lambdaSun += 2.0 * np.pi

    raSun, decSun = ecliptic2equatorial(lambdaSun, 0.0)

    return raSun, decSun






def skyToPlatformCoordinates(raStar, decStar, raPlatform, decPlatform, solarPanelOrientation):

    """From equatorial to platform coordinates.

    Convert the equatorial sky coordinates (alpha, delta) of a star 
    to cartesian platform coordinates (xSC, ySC, zSC).
    
    NOTE: Reference documents: PLATO-KUL-PL-TN-0001

    Parameters
    ----------
    raStar : float
        Right ascension of the star [rad]
    decStar : float   
        Declination of the star [rad]
    raPlatform : float
        Right ascension of the platform roll axis [rad]
    decPlatform : float
        Declination of the platform roll axis [rad]
    solarPanelOrientation : float
        Orientation of the solar panel [rad]
        This corresponds to (0, pi/2, pi, 3pi/2) for quarters (Q1,Q2,Q3,Q4)

    Return
    ------
    xSC, ySC, zSC : float  
        Normalized cartesian coordinates of the direction of the star
        in the spacecraft reference frame.
    """

    # Get the sky position of the Sun (ra, dec) [rad]

    raSun, decSun = sunSkyCoordinatesAwayfromPlatformPointing(raPlatform, decPlatform,
                                                              solarPanelOrientation)

    # Compute the equatorial cartesian coordinates of the unit vector
    # along the z-axis (= roll = pointing axis) of the platform.
    # The x-axis of the platform points to the highest point fof the
    # sunshield, which is pointing to the (average) sky position of the Sun.

    zSC = np.array([np.cos(decPlatform)*np.cos(raPlatform),
                    np.cos(decPlatform)*np.sin(raPlatform),
                    np.sin(decPlatform)])
    deltax = np.arctan(- np.cos(raPlatform-raSun) / np.tan(decPlatform))
    xSC = np.array([np.cos(deltax)*np.cos(raSun), np.cos(deltax)*np.sin(raSun), np.sin(deltax)])
    ySC = np.cross(zSC, xSC)

    # Compute the rotation matrix to convert cartesian coordinates in
    # the equatorial reference frame to cartesian coordinates in the spacecraft framework

    rotEQ2SC   = np.array([[xSC[0], xSC[1], xSC[2]], \
                           [ySC[0], ySC[1], ySC[2]], \
                           [zSC[0], zSC[1], zSC[2]]])

    # Compute the cartesian coordinates of the star in the equatorial reference frame

    starEQ = np.array([np.cos(decStar)*np.cos(raStar), np.cos(decStar)*np.sin(raStar), np.sin(decStar)])

    # Transform these coordinates to the corresponding ones in the focal plane reference frame

    starSC = np.dot(rotEQ2SC, starEQ)

    # That's it

    return starSC[0], starSC[1], starSC[2]





def platformToSkyCoordinates(xPLM, yPLM, zPLM, raPlatform, decPlatform, solarPanelOrientation):

    """From platoform pointing to equatorial sky pointing.

    Convert the cartesian platform coordinates (xPLM, yPLM, zPLM) of a point to equatorial
    sky coordinates (alpha, delta). The units of the platform coordinates are arbitrary,
    since the output of this function is only a sky 'direction'.

    NOTE: Reference documents: PLATO-KUL-PL-TN-0001

    Parameters
    ----------
    xPLM : float
        X-coordinate in the platoform reference frame [arbitrary]
    yPLM : float
        Y-coordinate in the platoform reference frame [same unit as xPLM]
    zPLM : float                     
        Z-coordinate in the platoform reference frame [same unit as xPLM]
    raPlatform : float
        Right ascension of the platform roll axis [rad]
    decPlatform : float            
        Declination of the platform roll axis [rad]
    solarPanelOrientation : float
        Solar panel orientation [rad]
        This corresponds to (0,pi/2,pi,3pi/2) for quarters (Q1,Q2,Q3,Q4)

    Return
    ------
    raStar : float
        Right ascension according to direction of the vector (xPLM, yPLM, zPLM) [rad]
    decStar : float
        Declination according to direction of the vector (xPLM, yPLM, zPLM) [rad]
    """

    vecPLM = np.array([xPLM, yPLM, zPLM])

    # Get the sky position of the Sun (ra, dec) [rad]

    raSun, decSun = sunSkyCoordinatesAwayfromPlatformPointing(raPlatform, decPlatform,
                                                              solarPanelOrientation)

    # Compute the equatorial cartesian coordinates of the unit vector
    # along the z-axis (= roll = pointing axis) of the platform.
    # The x-axis of the platform points to the highest point fof the sunshield,
    # which is pointing to the (average) sky position of the Sun.

    zPLM = np.array([np.cos(decPlatform)*np.cos(raPlatform),
                     np.cos(decPlatform)*np.sin(raPlatform),
                     np.sin(decPlatform)])
    deltax = np.arctan(- np.cos(raPlatform-raSun) / np.tan(decPlatform))
    xPLM = np.array([np.cos(deltax)*np.cos(raSun), np.cos(deltax)*np.sin(raSun), np.sin(deltax)])
    yPLM = np.cross(zPLM, xPLM)

    # Compute the rotation matrix to convert cartesian coordinates in the equatorial
    # reference frame to cartesian coordinates in the platform reference frame

    rotPLM2EQ = np.array([[xPLM[0], yPLM[0], zPLM[0]], \
                          [xPLM[1], yPLM[1], zPLM[1]], \
                          [xPLM[2], yPLM[2], zPLM[2]]])

    # Transform the unnormalized focal plane coordinates to the
    # corresponding ones in the equatorial reference frame

    vecEQ = np.dot(rotPLM2EQ, vecPLM)

    # Convert the cartesian equatorial coordinates to equatorial sky coordinates

    norm    = np.sqrt(vecEQ[0]*vecEQ[0] + vecEQ[1]*vecEQ[1] + vecEQ[2]*vecEQ[2])
    decStar = np.pi/2.0 - np.arccos(vecEQ[2]/norm);
    raStar  = np.arctan2(vecEQ[1], vecEQ[0]);

    # Ensure that the right ascension is positive

    if isinstance(raStar, np.ndarray):
        raStar[raStar < 0.0] += 2.*np.pi
    else:
        if (raStar < 0.0):
            raStar += 2.*np.pi

    # That's it!

    return raStar, decStar





def skyToFocalPlaneCoordinates(raStar, decStar,
                               raPlatform, decPlatform, solarPanelOrientation,
                               tiltAngle, azimuthAngle,
                               focalPlaneAngle, focalLength):

    """From equatorial to focal plane reference system.

    Convert the equatorial sky coordinates (alpha, delta) of a star to 
    undistorted normalized focal plane coordinates (xFP, yFP), assuming
    a spherical pinhole camera.

    Notes
    ----- 
    - The unit of the cartesian focal plane coordinates is the same as the one
      of focalLength. focalLength can be expressed in e.g. [mm] or in [pixels]. 
      If focalLength == 1.0, then the corresponding focal plane coordinates are
      called "normalized coordinates."
    - Reference documents: PLATO-KUL-PL-TN-0001 and PLATO-DLR-PL-TN-016

    Parameters
    ----------
    raStar : float
        Right ascension of the star [rad]
    decStar : float   
        Declination of the star [rad]
    raPlatform : float
        Right ascension of the platform roll axis [rad]
    decPlatform : float
        Declination of the platform roll axis [rad]
    solarPanelOrientation : float
        Orientation of the solar panel [rad]
        This corresponds to (0, pi/2, pi, 3pi/2) for quarters (Q1,Q2,Q3,Q4)
    tiltAngle : float
        Tilt (altitude) angle of the camera w.r.t. the platform z-axis [rad]
    azimuthAngle : float
        Azimuth angle of the camera on the platform [rad]
    focalPlaneAngle : float        
        Angle between the Y_CAM axis and the Y_FP axis: gamma_FP [rad]
    focalLength : float
        Focal length of the camera. Unit: see remarks.

    Return
    ------
    xFP : float
        Normalized x-coordinate of the project star in the FP reference frame.
    yFP : float
        Normalized y-coordinate of the project star in the FP reference frame.
    """

    # Get the sky position of the Sun (ra, dec) [rad]

    raSun, decSun = sunSkyCoordinatesAwayfromPlatformPointing(raPlatform,
                                                              decPlatform,
                                                              solarPanelOrientation)

    # Compute the equatorial cartesian coordinates of the unit vector
    # along the z-axis (= roll = pointing axis) of the platform.
    # The x-axis of the platform points to the highest point fof the sunshield,
    # which is pointing to the (average) sky position of the Sun.

    zPLM = np.array([np.cos(decPlatform)*np.cos(raPlatform),
                     np.cos(decPlatform)*np.sin(raPlatform),
                     np.sin(decPlatform)])
    deltax = np.arctan(-np.cos(raPlatform-raSun) / np.tan(decPlatform))
    xPLM   = np.array([np.cos(deltax)*np.cos(raSun), np.cos(deltax)*np.sin(raSun), np.sin(deltax)])
    yPLM   = np.cross(zPLM, xPLM)

    # Compute the rotation matrix to convert cartesian coordinates in the equatorial
    # reference frame to cartesian coordinates in the platform framework.

    rotEQ2PLM = np.array([[xPLM[0], xPLM[1], xPLM[2]],
                          [yPLM[0], yPLM[1], yPLM[2]],
                          [zPLM[0], zPLM[1], zPLM[2]]])

    # Compute the rotation matrix to convert cartesian coordinates in the platform
    # reference frame to cartesian coordinates in the camera reference frame

    rotAzimuth = np.array([[ np.cos(azimuthAngle), np.sin(azimuthAngle), 0],
                           [-np.sin(azimuthAngle), np.cos(azimuthAngle), 0],
                           [        0        ,          0,         1]])

    rotTilt = np.array([[np.cos(tiltAngle), 0, -np.sin(tiltAngle)],
                        [     0        , 1,        0       ],
                        [np.sin(tiltAngle), 0,  np.cos(tiltAngle)]])

    rotPLM2CAM = np.dot(rotAzimuth.T, np.dot(rotTilt, rotAzimuth))

    # Compute the rotation matrix to convert cartesian coordinates in the camera
    # reference frame to cartesian coordinates in the focal plane reference frame

    rotCAM2FP = np.array([[ np.cos(focalPlaneAngle), np.sin(focalPlaneAngle), 0],
                          [-np.sin(focalPlaneAngle), np.cos(focalPlaneAngle), 0],
                          [          0          ,           0         , 1]])

    # Combine all the rotation matrices

    rotEQ2FP = np.dot(rotCAM2FP, np.dot(rotPLM2CAM, rotEQ2PLM))

    # Compute the cartesian coordinates of the star in the equatorial reference frame

    starEQ = np.array([np.cos(decStar)*np.cos(raStar), np.cos(decStar)*np.sin(raStar), np.sin(decStar)])

    # Transform these coordinates to the corresponding ones in the focal plane reference frame:

    starFP = np.dot(rotEQ2FP, starEQ)

    # Convert the units to the one of focalLength (usually [mm]), and normalize the coordinates
    # to take into account the pinhole camera projection.

    xFPmm = -focalLength * starFP[0]/starFP[2]
    yFPmm = -focalLength * starFP[1]/starFP[2]

    # That's it

    return xFPmm, yFPmm





def focalPlaneToSkyCoordinates(xFP, yFP, raPlatform, decPlatform, solarPanelOrientation,
                               tiltAngle, azimuthAngle, focalPlaneAngle, focalLength):

    """From focal plane to equatorial reference system.

    Convert the undistorted normalized focal plane coordinates (xFP, yFP)
    of a star to equatorial sky coordinates (alpha, delta), assuming a 
    spherical pinhole camera.

    Notes
    -----
    The transformation assumes that the pinhole reverses the image.

    Parameters
    ----------
    xFP : float
        Undistorted normalized x-coordinate in the FP reference frame [unit as focalLength]
    yFP : float
        Undistorted normalized y-coordinate in the FP reference frame [unit as focalLength]
    raPlatform : float
        Right ascension of the platform roll axis [rad]
    decPlatform : float            
        Declination of the platform roll axis [rad]
    solarPanelOrientation : float
        Solar panel orientation [rad]
        This corresponds to (0,pi/2,pi,3pi/2) for quarters (Q1,Q2,Q3,Q4)
    tiltAngle : float
        Tilt (altitude) angle of the camera w.r.t. the platform z-axis [rad]
    azimuthAngle : float
        Azimuth angle of the camera on the platform [rad]
    focalPlaneAngle : float        
        Angle between the Y_CAM axis and the Y_FP axis: gamma_FP [rad]
    focalLength : float
        Focal length of the camera. Unit: see remarks.

    Return
    ------
    raStar : float
        Right ascension sky coordinate [rad]
    decStar : float
        Declination sky coordinate [rad]
    """

    # Get the sky position of the Sun (ra, dec) [rad]

    raSun, decSun = sunSkyCoordinatesAwayfromPlatformPointing(raPlatform,
                                                              decPlatform,
                                                              solarPanelOrientation)
    # Undo the reverse-image projection effect of the pinhole

    vecFP = np.array([-xFP/focalLength, -yFP/focalLength, 1.0], dtype=object)

    # Compute the rotation matrix to convert cartesian coordinates in the focal
    # plane reference frame to cartesian coordinates in the telescope reference frame

    rotFP2CAM = np.array([[np.cos(focalPlaneAngle), -np.sin(focalPlaneAngle), 0],
                          [np.sin(focalPlaneAngle),  np.cos(focalPlaneAngle), 0],
                          [          0         ,            0         , 1]])

    # Compute the rotation matrix to convert cartesian coordinates in the telescope
    # reference frame to cartesian coordinates in the spacecraft reference frame

    rotAzimuth = np.array([[np.cos(azimuthAngle), -np.sin(azimuthAngle), 0],
                           [np.sin(azimuthAngle),  np.cos(azimuthAngle), 0],
                           [        0        ,          0        , 1]])

    rotTilt = np.array([[ np.cos(tiltAngle), 0, np.sin(tiltAngle)],
                        [     0         , 1,        0      ],
                        [-np.sin(tiltAngle), 0, np.cos(tiltAngle)]])

    rotCAM2PLM = np.dot(rotAzimuth, np.dot(rotTilt, rotAzimuth.T))

    # Compute the equatorial cartesian coordinates of the unit vector along
    # the z-axis (= roll = pointing axis) of the platform.
    # The x-axis of the platform points to the highest point fof the sunshield,
    # which is pointing to the (average) sky position of the Sun.

    zPLM = np.array([np.cos(decPlatform)*np.cos(raPlatform),
                     np.cos(decPlatform)*np.sin(raPlatform),
                     np.sin(decPlatform)])
    deltax = np.arctan(- np.cos(raPlatform-raSun) / np.tan(decPlatform))
    xPLM = np.array([np.cos(deltax)*np.cos(raSun), np.cos(deltax)*np.sin(raSun), np.sin(deltax)])
    yPLM = np.cross(zPLM, xPLM)

    # Compute the rotation matrix to convert cartesian coordinates in the equatorial
    # reference frame to cartesian coordinates in the spacecraft reference frame

    rotPLM2EQ = np.array([[xPLM[0], yPLM[0], zPLM[0]], \
                          [xPLM[1], yPLM[1], zPLM[1]], \
                          [xPLM[2], yPLM[2], zPLM[2]]])

    # Combine all the rotation matrices

    rotFP2EQ = np.dot(rotPLM2EQ, np.dot(rotCAM2PLM, rotFP2CAM))

    # Transform the unnormalized focal plane coordinates to the corresponding
    # ones in the equatorial reference frame

    vecEQ = np.dot(rotFP2EQ, vecFP)

    # Convert the cartesian equatorial coordinates to equatorial sky coordinates

    norm    = np.sqrt(vecEQ[0]*vecEQ[0] + vecEQ[1]*vecEQ[1] + vecEQ[2]*vecEQ[2])
    decStar = np.pi/2.0 - np.arccos(vecEQ[2]/norm);
    raStar  = np.arctan2(vecEQ[1], vecEQ[0]);

    # Ensure that the right ascension is positive

    if isinstance(raStar, np.ndarray):
        raStar[raStar < 0.0] += 2.*np.pi
    else:
        if (raStar < 0.0):
            raStar += 2.*np.pi

    # That's it!

    return raStar, decStar





def pixelCoordinates2FocalPlaneAngles(xCCD, yCCD, ccdCode, pixelSize, focalLength):

    """From pixel to focal plane coordinates.

    Given the real-valued CCD pixel coordinates, compute the location angles
    of the star in the focal plane. These calculations are based on the custom
    CCD position and orientation angle (so not from file!).

    Parameters
    ----------
    xCCD : int
        Real-valued x-coordinate on the CCD (column) [pix]
    yCCD : int 
        Real-valued y-coordinates on the CCD (row) [pix]
    ccdCode : str
        For N-CAMs use: '1', '2', '3', '4'
        For F-CAMs use: '1F','2F','3F','4F'
    pixelSize : float
        Size of one square pixel [microns]
    focalLength : float
        Focal length of the camera [mm]

    Return
    ------
    angleFromOpticalAxis : float
        Angular distance from the optical axis [rad]
    azimuthFromXAxis : float 
        Azimuth angle from the X-axis of the reference CCD [rad]
    """

    ccdZeroPointX = CCD[ccdCode]["zeroPointXmm"]
    ccdZeroPointY = CCD[ccdCode]["zeroPointYmm"]
    ccdAngle      = CCD[ccdCode]["angle"]

    xmm, ymm = pixelToFocalPlaneCoordinates(xCCD, yCCD, pixelSize,
                                            ccdZeroPointX, ccdZeroPointY, ccdAngle)
    angleFromOpticalAxis = gnomonicRadialDistanceFromOpticalAxis(xmm,ymm,focalLength)
    azimuthFromXAxis = np.arctan2(ymm,xmm)

    return angleFromOpticalAxis, azimuthFromXAxis





def focalPlaneAngles2pixelCoordinates(angleFromOpticalAxis, azimuthFromXAxis,
                                      ccdCode, pixelSize, focalLength):

    """From focal plane angles to pixel coordinates.

    Given the location angles of the star in the focal plane, compute
    the real-valued CCD pixel coordinates. These calculations are based
    on the  custom CCD position and orientation angle (so not from file!).

    Parameters
    ----------
    angleFromOpticalAxis : float
        Angular distance from the optical axis [rad]
    azimuthFromXAxis : float 
        Azimuth angle from the X-axis of the reference CCD [rad]
    ccdCode : str
        For N-CAMs use: '1', '2', '3', '4'
        For F-CAMs use: '1F','2F','3F','4F'
    pixelSize : float
        Size of one square pixel [microns]
    focalLength : float
        Focal length of the camera [mm]

    Return
    ------
    xCCD : int
        Real-valued x-coordinate on the CCD (column) [pix]
    yCCD : int 
        Real-valued y-coordinates on the CCD (row) [pix]
    """

    ccdZeroPointX = CCD[ccdCode]["zeroPointXmm"]
    ccdZeroPointY = CCD[ccdCode]["zeroPointYmm"]
    ccdAngle      = CCD[ccdCode]["angle"]

    xFPmm, yFPmm = focalPlaneCoordinatesFromGnomonicRadialDistance(angleFromOpticalAxis,
                                                                   focalLength,
                                                                   inPlaneRotation=azimuthFromXAxis)
    xCCD, yCCD = focalPlaneToPixelCoordinates(xFPmm, yFPmm, pixelSize,
                                              ccdZeroPointX, ccdZeroPointY, ccdAngle)

    # That's it!
    
    return xCCD, yCCD





def telescopeToUndistortedFocalPlaneCoordinates(xCAM, yCAM, zCAM, focalLength, focalPlaneAngle):

    """From 'telescope' to undistorted focal plane coordinates.

    Given cartesian coordinates in the telescope reference frame, compute
    the undistorted coordinates in the focal plane reference frame.

    Parameters
    ----------
    xCAM : float
        X-coordinate in the camera reference frame [arbitrary unit]
    yCAM : float
        Y-coordinate in the camera reference frame [same unit as xPLM]
    zCAM : float
        Z-coordinate in the camera reference frame [same unit as xPLM]
    focalLength : float
        Focal length of the camera. Unit: see remarks.
    focalPlaneAngle : float        
        Angle between the Y_CAM axis and the Y_FP axis: gamma_FP [rad]

    Return
    ------
    xFP : float
        Undistorted normalized x-coordinate in the FP reference frame [unit as focalLength]
    yFP : float
        Undistorted normalized y-coordinate in the FP reference frame [unit as focalLength]
    """

    rotCAM2FP = np.array([[ np.cos(focalPlaneAngle), np.sin(focalPlaneAngle), 0],
                          [-np.sin(focalPlaneAngle), np.cos(focalPlaneAngle), 0],
                          [          0          ,           0         , 1]])

    vecCAM = np.array([xCAM, yCAM, zCAM])
    vecFP = np.dot(rotCAM2FP, vecCAM)

    # Convert the units to the one of focalLength (usually [mm]), and normalize the coordinates
    # to take into account the pinhole camera projection.

    xFPmm = -focalLength * vecFP[0]/vecFP[2]
    yFPmm = -focalLength * vecFP[1]/vecFP[2]

    return xFPmm, yFPmm





def undistortedFocalPlaneToTelescopeCoordinates(xFP, yFP, focalLength, focalPlaneAngle):

    """From undistorted focal plane to camera coordinates.
    
    Given cartesian coordinates in the focal plane reference frame, compute 
    the coordinates in the telescope reference frame.

    Notes
    -----
    Because of the projection degeneracy, the camera coordinates are
    computed on the unit sphere.

    Parameters
    ----------
    xFP : float
        Undistorted normalized x-coordinate in the FP reference frame [unit as focalLength]
    yFP : float
        Undistorted normalized y-coordinate in the FP reference frame [unit as focalLength]
    focalLength : float
        Focal length of the camera. Unit: see remarks.
    focalPlaneAngle : float        
        Angle between the Y_CAM axis and the Y_FP axis: gamma_FP [rad]

    Return
    ------
    xCAM : float
        X-coordinate in the camera reference frame [arbitrary]
    yCAM : float
        Y-coordinate in the camera reference frame [same unit as xCAM]
    zCAM : float                     
        Z-coordinate in the camera reference frame [same unit as xCAM]
    """

    vecFP = np.array([-xFP/focalLength, -yFP/focalLength, 1.0])

    # Compute the rotation matrix to convert cartesian coordinates in the focal plane
    # reference frame to cartesian coordinates in the telescope reference frame

    rotFP2CAM = np.array([[np.cos(focalPlaneAngle), -np.sin(focalPlaneAngle), 0],
                          [np.sin(focalPlaneAngle),  np.cos(focalPlaneAngle), 0],
                          [          0         ,            0         , 1]])

    vecCAM = np.dot(rotFP2CAM, vecFP)
    vecCAM /= norm(vecCAM)

    return vecCAM[0], vecCAM[1], vecCAM[2]





def undistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm, distortionCoefficients, focalLength):

    """From undistorted to distorted focal plane coordinates.

    Convert from undistorted to distorted normalized focal plane coordinates
    using the analytic distortion model.

    Notes
    -----
    - This is not the prefered method for detectors with mapped PSF, since these
      use a mapped distortion model. For such models use the function: 
      mappedUndistortedToDistortedFocalPlaneCoordinates
    - Example of distortion coefficients: 
      [0.32419, 0.0232909, 0.407979, 0.00022463, 0.000217599, 0.000381958, 0.000963902]

    Parameters
    ----------
    xFPmm : float 
        Undistorted focal plane x-coordinate [mm]
    yFPmm : float
        Undistorted focal plane y-coordinate [mm]
    distortionCoefficients : list, ndarray
        List of polynomial coefficients
    focalLength : float
        Focal length [mm]

    Return
    ------
    xFPdist : float
        Distorted focal plane x-coordinates [mm]
    yFPdist : float
        Distorted focal plane y-coordinates [mm]
    """

    radialCoefficients = [0, 0, 0,
                          distortionCoefficients[0], 0,
                          distortionCoefficients[1], 0,
                          distortionCoefficients[2]]
    distortionPolynomial = Polynomial(radialCoefficients)

    # Position angle on the focal plane [radians]
    
    angle = np.arctan2(yFPmm, xFPmm)

    # Undistorted radial distance [normalised pixels]
    
    rFP = np.sqrt(xFPmm**2 + yFPmm**2) / focalLength
    radialDistortion = distortionPolynomial(rFP) * focalLength
    tangentialDistortion = rFP**2 * (distortionCoefficients[5] * np.cos(angle) +
                                     distortionCoefficients[6] * np.sin(angle)) * focalLength

    xFPdist = (xFPmm + np.cos(angle) *
               (radialDistortion + tangentialDistortion) +
               distortionCoefficients[3] * rFP**2 * focalLength)
    yFPdist = (yFPmm + np.sin(angle) *
               (radialDistortion + tangentialDistortion) +
               distortionCoefficients[4] * rFP**2 * focalLength)

    # That's it!
    
    return xFPdist, yFPdist






def mappedUndistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm, pathToPsfFile):

    """From mapped undistorted to distorted focal plane coordinates.

    Convert from undistorted to distorted normalized focal plane
    coordinates using a mapped distortion model.

    Notes
    -----
    This is the prefered method for detectors with mapped PSF.
    For detectors that use analytic PSF use the function:
    undistortedToDistortedFocalPlaneCoordinates

    Parameters
    ----------
    xFPmm : float
        Undistorted x-coordinate in the FP reference frame [mm]
    yFPmm : float
        Undistorted y-coordinate in the FP reference frame [mm]
    pathToFoPsfFile : str
        Absolute path to the PSF file (HDF5) that the Coordinates map

    Return
    ------
    xFPdist : float
        Distorted focal plane x-coordinate [mm]
    yFPdist : float
        Distorted focal plane y-coordinate [mm]
    """

    # If radial distance is above 85 mm
    
    if (xFPmm**2 + yFPmm**2 > 85):
        return xFPmm, yFPmm

    # Check the path to the PSF file excists

    if not os.path.exists(pathToPsfFile):
        if os.path.exists(os.environ["PLATO_PROJECT_HOME"] + "/" + pathToPsfFile):
            pathToPsfFile = os.environ["PLATO_PROJECT_HOME"] + "/" + pathToPsfFile
        else:
            print(f"Error: {pathToPsfFile} is not a valid path name for mapped PSF")

    # We open the psf file where the coordinate transformation matrix should be in.
    # If the matrix isn't in the file, we raise an error.

    psfFile = h5py.File(pathToPsfFile, "r")
    if not "Coordinates map" in psfFile.keys():
        print("Error: No transformation map given in psf file, mapped distortion is not possible.")
        return
    else:
        coordMap = psfFile["Coordinates map"]

    # Calculate the distance of the undist coordinates wrt input coordinates

    undistorted = coordMap["Undistorted"]

    x = undistorted["x"]
    xUndis = np.zeros(x.shape, x.dtype)
    x.read_direct(xUndis)

    y = undistorted["y"]
    yUndis = np.zeros(y.shape, y.dtype)
    y.read_direct(yUndis)

    distorted = coordMap["Distorted"]

    x = distorted["x"]
    xDist = np.zeros(x.shape, x.dtype)
    x.read_direct(xDist)

    y = distorted["y"]
    yDist = np.zeros(y.shape, y.dtype)
    y.read_direct(yDist)
    
    distanceFromPointx = np.array([x - xFPmm for x in xUndis])
    distanceFromPointy = np.array([y - yFPmm for y in yUndis])
    aDistanceFromPoint  = np.array([ x**2 + y**2
        for x, y in zip(distanceFromPointx, distanceFromPointy)])

    # We should select the closest four undistorted points to the input point

    idx = np.arange(len(distanceFromPointx))
    idx_selected = np.empty(4, dtype=np.int16)

    idx_left  = idx[distanceFromPointx < 0]
    idx_right = idx[distanceFromPointx >= 0]

    leftDistanceFromPointy  = distanceFromPointy[distanceFromPointx < 0]
    rightDistanceFromPointy = distanceFromPointy[distanceFromPointx >= 0]

    left_bottom_idx = idx_left[leftDistanceFromPointy < 0]
    idx_closest_idx = np.argmin(aDistanceFromPoint[left_bottom_idx])
    idx_selected[0] = left_bottom_idx[idx_closest_idx]

    left_top_idx    = idx_left[leftDistanceFromPointy >=0]
    idx_closest_idx = np.argmin(aDistanceFromPoint[left_top_idx])
    idx_selected[1] = left_top_idx[idx_closest_idx]

    right_bottom_idx = idx_right[rightDistanceFromPointy < 0]
    idx_closest_idx  = np.argmin(aDistanceFromPoint[right_bottom_idx])
    idx_selected[2]  = right_bottom_idx[idx_closest_idx]

    right_top_idx = idx_right[rightDistanceFromPointy >= 0]
    idx_closest_idx = np.argmin(aDistanceFromPoint[right_top_idx])
    idx_selected[3] = right_top_idx[idx_closest_idx]

    for i in np.arange(2):
        if (yUndis[idx_selected[2*i]] > yUndis[idx_selected[2*i+1]]):
            dummy = idx_selected[2*i]
            idx_selected[2*i] = idx_selected[2*i+1]
            idx_selected[2*i+1] = dummy

    closestX = np.array([ xUndis[i] for i in idx_selected])
    closestY = np.array([ yUndis[i] for i in idx_selected])

    # We can write the points (xFPmm, yFPmm) as a linear combination of the
    # four closests points around this point.

    oPointIdx = [3, 2, 1, 0]

    constants = [abs( (closestX[oPointIdx[i]] - xFPmm) *
                      (closestY[oPointIdx[i]] - yFPmm))
                 for i in np.arange(4)]

    constants = [ constant / sum(constants) for constant in constants]

    closestXdist = np.array([ xDist[i] for i in idx_selected])
    closestYdist = np.array([ yDist[i] for i in idx_selected])

    xFPdist = sum([constants[i]*closestXdist[i] for i in np.arange(4)])
    yFPdist = sum([constants[i]*closestYdist[i] for i in np.arange(4)])

    # That's it!
    
    return xFPdist, yFPdist





def distortedToUndistortedFocalPlaneCoordinates(xFPdist, yFPdist,
                                                inverseDistortionCoefficients,
                                                focalLength):

    """From distorted to undistorted focal plane coordinates.

    Convert from distorted to undistorted normalized focal plane coordinates
    using the analytic distortion model.

    Notes
    -----
    - This is not the prefered method for detectors with mapped PSF, since
      these use a mapped distortion model. For such models use the function: 
      mappedDistortedToUndistortedFocalPlaneCoordinates.
    - Example of inverse distortion coefficients: 
      [-0.323487, 0.268344, -0.435473, -0.00019304, -0.000176961, -0.000321713, -0.000827654]

    Parameters
    ----------
    xFPdist : float
        Distorted focal plane x-coordinates [mm]
    yFPdist : float
        Distorted focal plane y-coordinates [mm]
    inverseDistortionCoefficients : list
        List of polynomial coefficients
    focalLength : float
        Focal length [mm]

    Return
    ------
    xFPmm : float
        Distorted x-coordinate in the FP reference frame [mm]
    yFPmm : float
        Distorted y-coordinate in the FP reference frame [mm]
    """

    inverseCoefficientsRadial = [0, 0, 0,
                                 inverseDistortionCoefficients[0], 0,
                                 inverseDistortionCoefficients[1], 0,
                                 inverseDistortionCoefficients[2]]
    inverseDistortionPolynomialRadial = Polynomial(inverseCoefficientsRadial)

    # Position angle on the focal plane [radians]
    
    angle = np.arctan2(yFPdist, xFPdist)

    # Distorted radial distance [normalised pixels]
    
    rFP = np.sqrt(xFPdist**2 + yFPdist**2) / focalLength

    # Distortion [mm] -> negative!
    
    radialDistortion     = inverseDistortionPolynomialRadial(rFP) * focalLength
    tangentialDistortion = (rFP**2 * (inverseDistortionCoefficients[5] * np.cos(angle) +
                                      inverseDistortionCoefficients[6] * np.sin(angle)) * focalLength)

    xFPmm = (xFPdist + np.cos(angle) * (radialDistortion + tangentialDistortion) +
             inverseDistortionCoefficients[3] * rFP**2 * focalLength)
    yFPmm = (yFPdist + np.sin(angle) * (radialDistortion + tangentialDistortion) +
             inverseDistortionCoefficients[4] * rFP**2 * focalLength)

    # That's it!
    
    return xFPmm, yFPmm





def mappedDistortedToUndistortedFocalPlaneCoordinates(xFPdist, yFPdist, pathToPsfFile):

    """From mapped distorted to undistorted focal plane coordinates.

    Convert from distorted to undistorted normalized focal plane coordinates
    using the mapped distortion model.

    Notes
    -----
    This is the prefered method for detectors with mapped PSF, since these use a
    mapped distortion model. For detectors that use analytic PSF use the function:
    distortedToUndistortedFocalPlaneCoordinates.

    Parameters
    ----------
    xFPdist : float
        Distorted focal plane x-coordinate [mm]
    yFPdist : float
        Distorted focal plane y-coordinate [mm]
    pathToFoPsfFile : hdf5 file
        Absolute path to the PSF file that the Coordinates map

    Return
    ------
    xFPmm : float
        Distorted x-coordinate in the FP reference frame [mm]
    yFPmm : float
        Distorted y-coordinate in the FP reference frame [mm]
    """

    delta  = 100.
    length = 80.
    x0 = 0
    y0 = 0
    i  = 0

    while((delta > .001) and (i<160)):
        xDist, yDist = mappedUndistortedToDistortedFocalPlaneCoordinates(x0,
                                                        y0, pathToPsfFile)

        length = 3*length / 5

        if (xFPdist > xDist):
            if ((x0 + length) > 85):
                x0 = 85
            else:
                x0 = x0 + length
        elif (xFPdist < xDist):
            if ((x0 - length) < -85):
                x0 = -85
            else:
                x0 = x0 - length

        if (yFPdist > yDist):
            if ((y0 + length) > 85):
                y0 = 85
            else:
                y0 = y0 + length
        elif (yFPdist < yDist):
            if ((y0 - length) < -85):
                y0 = -85
            else:
                y0 = y0 - length


        delta = abs(xFPdist - xDist) + abs(yFPdist - yDist)

        i += 1

    return x0, y0





def pixelToFocalPlaneCoordinates(xCCDpixel, yCCDpixel, pixelSize,
                                 ccdZeroPointX, ccdZeroPointY, CCDangle):

    """From pixel to focal plane coordinates.

    Given the (real-valued) pixel coordinates of the star on the CCD, 
    compute the (x,y) coordinates in the FP reference system.

    Parameters
    ----------
    xCCD : int
        X-coordinate of star on the CCD (column) [pix]
    yCCD : int 
        Y-coordinate of star on the CCD (row) [pix]
    pixelSize : float
        Size of one square pixel [microns]
    ccdZeroPointX : float
        X-coordinate of the CCD (0,0) point in the FP reference system [mm]
    ccdZeroPointY : float
        Y-coordinate of the CCD (0,0) point in the FP reference system [mm]
    CCDangle : float
        CCD orientation angle in the FP reference frame  [rad]

    Return
    ------
    xFP : float
        Column pixel coordinate of the star (real-valued) [mm]
    yFP : float
        Row pixel coordinate of the star (real-valued) [mm]
    """

    # Convert the pixel coordinates into [mm] coordinates

    xCCDmm = xCCDpixel * pixelSize / 1000.0
    yCCDmm = yCCDpixel * pixelSize / 1000.0

    # Convert the CCD coordinates into FP coordinates [mm]

    xFP = (xCCDmm - ccdZeroPointX) * np.cos(CCDangle) - (yCCDmm - ccdZeroPointY) * np.sin(CCDangle)
    yFP = (xCCDmm - ccdZeroPointX) * np.sin(CCDangle) + (yCCDmm - ccdZeroPointY) * np.cos(CCDangle)

    # That's it

    return xFP, yFP





def focalPlaneToPixelCoordinates(xFP, yFP, pixelSize, ccdZeroPointX, ccdZeroPointY, CCDangle):

    """From focal plaen to pixel coordinate.

    Compute the (real-valued) pixel coordinates of the star on the CCD, given the (x,y)
    coordinates in the FP reference system. The cartesian focal plane coordinates are
    supposed to be normalized (i.e. with the pinhole projection taken into account).

    Parameters
    ----------
    xFP : float
        Normalized x-coordinate of star in the FP reference frame [mm]
    yFP : float
        Normalized y-coordinate of star in the FP reference frame [mm]
    pixelSize : float
        Size of one square pixel [microns]
    ccdZeroPointX : float
        X-coordinate of the CCD (0,0) point in the FP reference system [mm]
    ccdZeroPointY : float
        Y-coordinate of the CCD (0,0) point in the FP reference system [mm]
    CCDangle : float
        CCD orientation angle in the FP reference frame  [rad]

    Return
    ------
    xCCDpixel : float
        Column pixel coordinate of the star (real-valued)
    yCCDpixel : float
        Rrow pixel coordinate of the star (real-valued)
    """

    # Convert the FP coordinates into CCD coordinates [mm]

    xCCDmm = ccdZeroPointX + xFP * np.cos(CCDangle) + yFP * np.sin(CCDangle)
    yCCDmm = ccdZeroPointY - xFP * np.sin(CCDangle) + yFP * np.cos(CCDangle)

    # Convert the [mm] coordinates into pixel coordinates

    xCCDpixel = xCCDmm / pixelSize * 1000.0
    yCCDpixel = yCCDmm / pixelSize * 1000.0

    # That's it

    return xCCDpixel, yCCDpixel





def gnomonicRadialDistanceFromOpticalAxis(xFP, yFP, focalLength):

    """Gnomonic radial distance away from optical axis.

    Calculate the gnomonic radial distance with respect to the 
    optical axis in the focal plane.

    Parameters
    ----------
    xFP : float
        Focal plane x-coordinate [mm]
    yFP : float
        Focal plane y-coordinate [mm]
    focalLength : float
        focal length of the camera [mm]

    Return
    ------
    angularDistance : float
        Rhe angular distance of the star w.r.t. the optical axis [rad]
    """

    tanx = xFP / focalLength
    tany = yFP / focalLength

    angularDistance = np.arccos(1.0/np.sqrt(1.0 + tanx*tanx + tany*tany));

    # Take care that the angle is between [0, 2*PI]

    if angularDistance < 0.0:
        angularDistance += 2.0 * np.pi
    elif angularDistance > 2.0 * np.pi:
        angularDistance -= 2.0 * np.pi

    # That's it!

    return angularDistance;





def focalPlaneCoordinatesFromGnomonicRadialDistance(angularDistance, focalLength,
                                                    inPlaneRotation=0.):

    """From focal plane coordinates to gnomonic radial distance.

    Calculate the xFP,yFP focal plane coordinates from the gnomonic
    radial distance with respect to the optical axis in the focal plane

    Parameters
    ----------
    angularDistance : float
        Rhe angular distance of the star w.r.t. the optical axis [rad]
    focalLength : float
        focal length of the camera [mm]
    inPlaneRotation : float
        Angle from the xFP axis to the target (default=0) [rad]

    Return
    ------
    xFP : float
        Focal plane x-coordinate [mm]
    yFP : float
        Focal plane y-coordinate [mm]
    """

    D = focalLength * np.tan(angularDistance)

    xFP = D * np.cos(inPlaneRotation)
    yFP = D * np.sin(inPlaneRotation)

    # That's it!
    
    return xFP,yFP





def computeCCDcornersInFocalPlane(ccdCode, pixelSize):

    """Compute CCD corners in focal plane.

    Get the (x,y) coordinates of each of the 4 corners of the exposed part
    of the CCD in the FP' reference system.  These calculations are based
    on the custom CCD position and orientation angle (so not from file!).

    Parameters
    ----------
    ccdCode : str
        For N-CAMs use: '1', '2', '3', '4'
        For F-CAMs use: '1F','2F','3F','4F'
    pixelSize : float
        Size of one square pixel [microns]

    Return
    ------
    cornersXmm : float
        X-coordinates of each of the corners in the FP' reference system [mm]
    cornersYmm : float
        Y-coordinates of each of the corners in the FP' reference system [mm]
    """

    # Get the pixel coordinates of the 4 corners of the exposed part of the CCD
    # Note that the x-direction corresponds to the CCD columns,
    # and the y-direction to the CCD rows.

    Nrows = CCD[ccdCode]["Nrows"]
    Ncols = CCD[ccdCode]["Ncols"]
    firstRow = CCD[ccdCode]["firstRow"]

    cornersXpix = array([0.0, Ncols, Ncols, 0.0])
    cornersYpix = array([firstRow, firstRow, Nrows, Nrows])

    # Convert to the x,y coordinates in the FP' reference frame

    zeroPointXmm = CCD[ccdCode]["zeroPointXmm"]
    zeroPointYmm = CCD[ccdCode]["zeroPointYmm"]
    ccdAngle     = CCD[ccdCode]["angle"]

    cornersXmm, cornersYmm = pixelToFocalPlaneCoordinates(cornersXpix, cornersYpix, pixelSize,
                                                          zeroPointXmm, zeroPointYmm, ccdAngle)

    # That's it

    return cornersXmm, cornersYmm






def getCCDandPixelCoordinates(raStar, decStar, raPlatform, decPlatform, solarPanelOrientation,
                              tiltAngle, azimuthAngle, focalPlaneAngle, focalLength, pixelSize,
                              includeFieldDistortion, normal, mappedDistortion=False,
                              distortionCoefficients=None, pathToPsfFile=None):

    """Get the CCD and pixel coordinates.

    Given the equatorial coordinates of a star, find out on which CCD
    it falls ('1', '2', '3', '4') and compute the pixel coordinates of
    the star on this CCD. If the star doesn't fall on any of the CCDs
    then (None, None, None) is given as output. These calculations use
    the custom CCD positions and orientation angle (so not from file!).

    Parameters
    ----------
    raStar : float
        Right ascension of the star [rad]
    decStar : float   
        Declination of the star [rad]
    raPlatform : float
        Right ascension of the platform roll axis [rad]
    decPlatform : float
        Declination of the platform roll axis [rad]
    solarPanelOrientation : float
        Orientation of the solar panel [rad]
        This corresponds to (0, pi/2, pi, 3pi/2) for quarters (Q1, Q2, Q3, Q4)
    tiltAngle : float
        Tilt (altitude) angle of the camera w.r.t. the platform z-axis [rad]
    azimuthAngle : float
        Azimuth angle of the camera on the platform [rad]
    focalPlaneAngle : float        
        Angle between the Y_CAM axis and the Y_FP axis: gamma_FP [rad]
    focalLength : float
        Focal length of the camera. Unit: see remarks.
    pixelSize : float
        Size of one square pixel [microns]
    includeFieldDistortion : bool
        True to include field distortion in coordinate transformations, false otherwise
    normal : bool
        True for the normal camera configuration, False for the fast cameras
    mappedDistortion : bool
        True if we want mapped distortion (mapped from file psf) False if we have analytic psfs
    distortionCoefficients : list
        Coefficients of the polynomial describing the distortion for anlytic psf
    pathToPsfFile : str
        Path to the PSF file (HDF5) for mapped PSFs to calculate mapped distortion

    Return
    ------
    ccdCode : str
        For N-CAMs use: '1', '2', '3', '4'
        For F-CAMs use: '1F','2F','3F','4F'
    xCCD : int
        X-coordinate (column) of star on the CCD [pix]
        If not on CCD: None
    yCCD : int 
        Y-coordinate (row) of star on the CCD (row) [pix]
        If not on CCD: None
    """

    # Make sure that for the respective field distortion the proper information is given.

    if (includeFieldDistortion or includeFieldDistortion == "yes"):
        if (mappedDistortion and pathToPsfFile is None):
            print("Error: If mapped field distortion should be taken into account, " +
                  "a path to the psf file should be given")
            return
        elif ( (not mappedDistortion) and distortionCoefficients is None):
            print("Error: If analytic field distortion should be taken into account, " +
                  "the distortionCoefficients should be given")
            return

    # Select the proper CCD codes depending on N-CAMs or F-CAMs

    if normal == True:
        ccdCodes = ['1', '2', '3', '4']
    else:
        ccdCodes = ['1F', '2F', '3F', '4F']

    # Compute the (x,y) coordinates in the FP reference system [mm]

    xFPmm, yFPmm = skyToFocalPlaneCoordinates(raStar, decStar,
                                              raPlatform, decPlatform, solarPanelOrientation,
                                              tiltAngle, azimuthAngle,
                                              focalPlaneAngle, focalLength)

    if (includeFieldDistortion == True) or (includeFieldDistortion == "yes"):
        if mappedDistortion:
            xFPmm, yFPmm = mappedUndistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm,
                                                                             pathToPsfFile)
        else:
            xFPmm, yFPmm = undistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm,
                                                                       distortionCoefficients,
                                                                       focalLength)

    # Find out if this falls on a CCD, and if yes which one.
    # Our approach: try each of the CCDs. Not elegant, but robust!

    for ccdCode in ccdCodes:

        # Compute the position of the star in pixel coordinates, for the current CCD,
        # disregarding the physical extend of the CCD

        zeroPointXmm = CCD[ccdCode]["zeroPointXmm"]
        zeroPointYmm = CCD[ccdCode]["zeroPointYmm"]
        ccdAngle     = CCD[ccdCode]["angle"]

        xCCDpix, yCCDpix = focalPlaneToPixelCoordinates(xFPmm, yFPmm, pixelSize,
                                                        zeroPointXmm, zeroPointYmm, ccdAngle)

        # Check if the star falls on the exposed area of the CCD. If not: go to next CCD

        Nrows = CCD[ccdCode]["Nrows"]
        Ncols = CCD[ccdCode]["Ncols"]
        firstRow = CCD[ccdCode]["firstRow"]

        if (xCCDpix < 0)      or (yCCDpix < firstRow): continue
        if (xCCDpix >= Ncols) or (yCCDpix >= Nrows):   continue

        # If we arrive here, we found a CCD on which the star is located

        return ccdCode, xCCDpix, yCCDpix

    # If we arrive here, the star does not fall on any CCD

    return None, None, None





def platformToTelescopePointingCoordinates(raPlatform, decPlatform, raSun, decSun,
                                           azimuthAngle, tiltAngle):

    """From platform to camera pointing coordinates.

    Given the platform pointing coordinates (i.e. the sky coordinates of the jitter axis)
    and the orientation of the telescope on the platform, compute the sky coordinates of
    the optical axis of the telescope. See also: PLATO-KUL-PL-TN-001.

    Parameters
    ----------
    raPlatform : float
        Right ascension of the pointing of the Platform [rad]
    decPlatform : float 
        Declination of the pointing of the Platform [rad]
    raSun : float
        Right ascension of the Sun (shield) at the given time [rad]
    decSun : float
        Declination of the Sun (shield) at the given time [rad]
    azimuthAngle : float
        Azimuth angle of the camera on the platform [rad]
    tiltAngle : float
        Tilt (altitude) angle between platform and camera pointing axes [rad]

    Return
    ------
    raTelescope : float
        Right ascension of the optical axis of the telescope [rad]
    decTelescope : float
        Declination of the optical axis of the telescope [rad]
    """

    # Compute the rotation matrix to convert cartesian coordinates in the telescope
    # reference frame to cartesian coordinates in the spacecraft reference frame

    rotAzimuth = np.array([[np.cos(azimuthAngle), -np.sin(azimuthAngle), 0],
                           [np.sin(azimuthAngle),  np.cos(azimuthAngle), 0],
                           [        0        ,          0,         1]])

    rotTilt = np.array([[ np.cos(tiltAngle), 0, np.sin(tiltAngle)], \
                        [     0        , 1,        0       ], \
                        [-np.sin(tiltAngle), 0, np.cos(tiltAngle)]])

    rotCAM2PLM = np.dot(rotAzimuth, np.dot(rotTilt, rotAzimuth.T))

    # Compute the equatorial cartesian coordinates of the unit vector
    # along the z-axis (= roll = pointing axis) of the platform.
    # The x-axis of the platform points to the highest point fof the sunshield,
    # which is pointing to the (average) sky position of the Sun.

    zPLM = np.array([np.cos(decPlatform)*np.cos(raPlatform),
                     np.cos(decPlatform)*np.sin(raPlatform),
                     np.sin(decPlatform)])
    deltax = np.arctan(- np.cos(raPlatform-raSun) / np.tan(decPlatform))
    xPLM = np.array([np.cos(deltax)*np.cos(raSun), np.cos(deltax)*np.sin(raSun), np.sin(deltax)])
    yPLM = np.cross(zPLM, xPLM)

    # Compute the rotation matrix to convert cartesian coordinates in the equatorial
    # reference frame to cartesian coordinates in the spacecraft reference frame

    rotPLM2EQ = np.array([[xPLM[0], yPLM[0], zPLM[0]],
                          [xPLM[1], yPLM[1], zPLM[1]],
                          [xPLM[2], yPLM[2], zPLM[2]]])

    # Combine all the rotation matrices

    rotCAM2EQ = np.dot(rotPLM2EQ, rotCAM2PLM)

    # In the telescope reference frame, the optical axis is simply the z-axis = (0,0,1)

    vecCAM = np.array([0.0, 0.0, 1.0])

    # Get the equatorial coordinates of this optical axis vector.

    vecEQ = np.dot(rotCAM2EQ, vecCAM)

    # Convert the cartesian equatorial coordinates to equatorial sky coordinates

    norm           = np.sqrt(vecEQ[0]*vecEQ[0] + vecEQ[1]*vecEQ[1] + vecEQ[2]*vecEQ[2])
    decOpticalAxis = np.pi/2.0 - np.arccos(vecEQ[2]/norm)
    raOpticalAxis  = np.arctan2(vecEQ[1], vecEQ[0])

    # Ensure that the right ascension is positive

    if isinstance(raOpticalAxis, np.ndarray):
        raOpticalAxis[raOpticalAxis < 0.0] += 2.*np.pi
    else:
        if (raOpticalAxis < 0.0):
            raOpticalAxis += 2.*np.pi

    # That's it

    return (raOpticalAxis, decOpticalAxis)





def getCameraGroupCoordinates(raPlatform, decPlatform, solarPanelOrientation=0):

    """Calculate the ICRS coordinates of each camera group.

    Parameters
    ----------
    raPlatform : float
        Right ascension of the pointing of the Platform [rad]
    decPlatform : float 
        Declination of the pointing of the Platform [rad]
    solarPanelOrientation : float
        Orientation of the solar panel [rad]
        This corresponds to (0, pi/2, pi, 3pi/2) for quarters (Q1,Q2,Q3,Q4)

    Return
    ------
    raGroups : float
        Right ascension of each camera group [rad]
    decGroups : float
        Declination of each camera group [rad]
    """

    # Find coordinates of the Sun [rad]

    raSun, decSun = sunSkyCoordinatesAwayfromPlatformPointing(np.deg2rad(raPlatform),
                                                              np.deg2rad(decPlatform),
                                                              np.deg2rad(solarPanelOrientation))

    # Relative azimuth and tilt angles of each camera group w.r.t platform

    azimuthAngles = [45.0, 135.0, 225.0, 315.0]
    tiltAngles    = [9.2, 9.2, 9.2, 9.2]

    # Fetch RA and Dec for each camera group

    raGroups  = []
    decGroups = []

    for group in range(4):
        ra, dec = platformToTelescopePointingCoordinates(np.deg2rad(raPlatform),
                                                         np.deg2rad(decPlatform),
                                                         raSun, decSun,
                                                         np.deg2rad(azimuthAngles[group]),
                                                         np.deg2rad(tiltAngles[group]))
        raGroups.append(np.rad2deg(ra))
        decGroups.append(np.rad2deg(dec))

    # Finito!

    return raGroups, decGroups





def calculateSubfieldAroundCoordinates(subfieldSizeX, subfieldSizeY, raStar, decStar,
                                       raPlatform, decPlatform, solarPanelOrientation,
                                       tiltTelescope, azimuthTelescope,
                                       focalPlaneAngle, focalLength, pixelSize,
                                       includeFieldDistortion, normal,
                                       mappedDistortion=False, distortionCoefficients=None,
                                       pathToPsfFile=None ):

    """Calculate location of subfield around equatorial coordinates.

    Calculates the location of the subfield such that the star with
    coordinates (raStar, decStar) is centered in the subfield. The
    function also checks if there is enough space around the pixel
    (feature/feature-447-CCD).

    Notes
    -----
    - This function is used by setSubfieldAroundCoordinates() and usually
      does not need to be called by the user.
    - If the coordinates do not fall on any CCD, an error message is shown.
    - If the star is too close to the edge for the given subfield size, and
      error message is shown, followed by an exit(1)

    Parameters
    ----------
    subfieldSizeX : int
       Full width (# of columns) of the subfield [pix]
    subfieldSizeY : int
       Full height (#of rows) of the subfield [pix]
    raStar : float
        Right ascension of the star [rad]
    decStar : float   
        Declination of the star [rad]
    raPlatform : float
        Right ascension of the platform roll axis [rad]
    decPlatform : float
        Declination of the platform roll axis [rad]
    solarPanelOrientation : float
        Orientation of the solar panel [rad]
        This corresponds to (0, pi/2, pi, 3pi/2) for quarters (Q1, Q2, Q3, Q4)
    tiltAngle : float
        Tilt (altitude) angle of the camera w.r.t. the platform z-axis [rad]
    azimuthAngle : float
        Azimuth angle of the camera on the platform [rad]
    focalPlaneAngle : float        
        Angle between the Y_CAM axis and the Y_FP axis: gamma_FP [rad]
    focalLength : float
        Focal length of the camera. Unit: see remarks.
    pixelSize : float
        Size of one square pixel [microns]
    includeFieldDistortion : bool
        True to include field distortion in coordinate transformations, false otherwise
    normal : bool
        True for the normal camera configuration, False for the fast cameras
    mappedDistortion : bool
        True if we want mapped distortion (mapped from file psf) False if we have analytic psfs
    distortionCoefficients : list
        Coefficients of the polynomial describing the distortion for anlytic psf
    pathToPsfFile : str
        Path to the PSF file (HDF5) for mapped PSFs to calculate mapped distortion

    Return
    ------
    ccdCode : str
        For N-CAMs use: '1', '2', '3', '4'
        For F-CAMs use: '1F','2F','3F','4F'
    xCCDpix : int
        X-coordinate (column) of star on the CCD [pix]
        If not on CCD: None
    yCCDpix : int 
        Y-coordinate (row) of star on the CCD [pix]
        If not on CCD: None
    """

    # Find out that we have been given the correct distortion input parameters.
    #If this is not the case raise error and return.
    
    if (includeFieldDistortion or includeFieldDistortion == "yes"):
        if (mappedDistortion and pathToPsfFile is None):
            print("Error: If mapped field distortion should be taken into account, " +
                  "a path to the psf file should be given")
            return
        elif ( (not mappedDistortion) and distortionCoefficients is None):
            print("Error: If analytic field distortion should be taken into account, " +
                  "the distortionCoefficients should be given")
            return

    # Find out on which CCD the star falls, and the corresponding pixel coordinates

    ccdCode, xCCDpix, yCCDpix = getCCDandPixelCoordinates(raStar, decStar,
                                                          raPlatform, decPlatform,
                                                          solarPanelOrientation,
                                                          tiltTelescope, azimuthTelescope,
                                                          focalPlaneAngle, focalLength,
                                                          pixelSize, includeFieldDistortion, normal,
                                                          mappedDistortion, distortionCoefficients,
                                                          pathToPsfFile)

    # If the CCD code is None, the star does not fall on any ccd -> error

    if ccdCode == None:
        return None, None, None

    # If the star does fall on a CCD, check if it's not too close to the edge for the subfield to
    # be completely on the CCD.

    xCCDpix = int(xCCDpix)               # integer values
    yCCDpix = int(yCCDpix)
    firstRow = CCD[ccdCode]["firstRow"]  # different from nominal than for fast cams
    Ncols = CCD[ccdCode]["Ncols"]
    Nrows = CCD[ccdCode]["Nrows"]

    if (xCCDpix - subfieldSizeX/2 < 0           or
        xCCDpix + subfieldSizeX/2 - 1 > Ncols-1 or
        yCCDpix - subfieldSizeY/2 < firstRow    or
        yCCDpix + subfieldSizeY/2 - 1 > Nrows-1):
        return None, None, None

    # That's it!

    return ccdCode, xCCDpix, yCCDpix






def skyToPixelCoordinates(sim, raStar, decStar, normal):

    """From equatorial to pixel coordinates.

    Convert sky coordinates to pixel coordinates. These calculations are based
    on the custom CCD position and orientation angle (so not from file!).

    Notes
    -----
    It is assumed that the configuration parameters in the sim object contains
    a correct (ra, dec) of the platform, a correct (azimuth, tilt) of the telescope,
    a valid value for the focal length, the plate scale, the pixel size, and that
    the switch to include distortion or not is set correctly.

    Parameters
    ----------
    sim : class instance
        Instance of Simulation class (simulation.py)
    raStar : float
        Right ascension of the star [rad]
    decStar : float   
        Declination of the star [rad]
    subfieldSizeX : int
       Full width (# of columns) of the subfield [pix]
    subfieldSizeY : int
       Full height (# of rows) of the subfield [pix]
    normal : bool
       True for the normal camera configuration, False for the fast cameras

    Return
    ------
    ccdCode : str
        For N-CAMs use: '1', '2', '3', '4'
        For F-CAMs use: '1F','2F','3F','4F'
    xCCDpix : int
        X-coordinate (column) of star on the CCD [pix]
        If not on CCD: None
    yCCDpix : int 
        Y-coordinate (row) of star on the CCD [pix]
        If not on CCD: None
    """

    # Resolve which distortion model is used (if any)
    
    if (sim["PSF/Model"] == "MappedFromFile"):
        includeFieldDistortion = True
        distortionCoefficients = None
        pathToPsfFile          = sim["PSF/MappedFromFile/Filename"]
        mappedDistortion       = True
    elif (sim["Camera/IncludeFieldDistortion"] == "yes"  or
          sim["Camera/IncludeFieldDistortion"] == "1"    or
          sim["Camera/IncludeFieldDistortion"]):
        distortionCoefficients = sim["Camera/FieldDistortion/ConstantCoefficients"]
        pathToPsfFile          = None
        mappedDistortion       = False
        includeFieldDistortion = True
    else:
        includeFieldDistortion = False
        distortionCoefficients = None
        pathToPsfFile          = None
        mappedDistortion       = False

    # Fetch parameters and convert
    
    pixelSize             = float(sim["CCD/PixelSize"])
    focalLength           = float(sim["Camera/FocalLength/ConstantValue"]) * 1000.0  # [m] -> [mm]
    raPlatform            = np.deg2rad(float(sim["ObservingParameters/RApointing"]))
    decPlatform           = np.deg2rad(float(sim["ObservingParameters/DecPointing"]))
    solarPanelOrientation = np.deg2rad(float(sim["Platform/SolarPanelOrientation"]))
    focalPlaneAngle       = np.deg2rad(float(sim["Camera/FocalPlaneOrientation/ConstantValue"]))
    azimuthTelescope      = np.deg2rad(float(sim["Telescope/AzimuthAngle"]))
    tiltTelescope         = np.deg2rad(float(sim["Telescope/TiltAngle"]))

    # Get the pixel coordinates on the CCD

    ccdCode, xCCDpixel, yCCDpixel = getCCDandPixelCoordinates(raStar, decStar,
                                                              raPlatform, decPlatform,
                                                              solarPanelOrientation,
                                                              tiltTelescope, azimuthTelescope,
                                                              focalPlaneAngle, focalLength,
                                                              pixelSize, includeFieldDistortion,
                                                              normal, mappedDistortion,
                                                              distortionCoefficients, pathToPsfFile)

    # That's it!
    
    return ccdCode, xCCDpixel, yCCDpixel






def pixelToSkyCoordinates(sim, ccdCode, xCCDpix, yCCDpix):

    """From pixel to equatorial coordinates.

    Convert pixel coordinates to equatorial sky coordinates (ra,dec).
    These calculations are based on the custom CCD position and 
    orientation angle (so not from file!).

    Notes
    -----
    It is assumed that the configuration parameters in the sim object 
    contains a correct (ra, dec) of the platform, a correct (azimuth, tilt)
    of the telescope, a valid value for the focal length, the pixel size,
    and that the switch to include distortion or not is set correctly.

    Parameters
    ----------
    sim : class instance
        Instance of Simulation class (simulation.py)
    ccdCode : str
        For N-CAMs use: '1', '2', '3', '4'
        For F-CAMs use: '1F','2F','3F','4F'
    xCCDpix : int
        X-coordinate (column) of star on the CCD [pix]
        If not on CCD: None
    yCCDpix : int 
        Y-coordinate (row) of star on the CCD [pix]
        If not on CCD: None

    Return
    ------
    raStar : float
        Right ascension of the star [rad]
    decStar : float   
        Declination of the star [rad]
    """

    # Resolve which distortion model is used (if any)
    
    if (sim["PSF/Model"] == "MappedFromFile"):
        includeFieldDistortion = True
        inverseDistortionCoefficients = None
        pathToPsfFile          = sim["PSF/MappedFromFile/Filename"]
        mappedDistortion       = True

    elif (sim["Camera/IncludeFieldDistortion"] == "yes" or
          sim["Camera/IncludeFieldDistortion"] == "1"   or
          sim["Camera/IncludeFieldDistortion"] == "True"):
        inverseDistortionCoefficients = sim["Camera/FieldDistortion/ConstantInverseCoefficients"]
        pathToPsfFile          = None
        mappedDistortion       = False
        includeFieldDistortion = True

    else:
        includeFieldDistortion        = False
        pathToPsfFile                 = None
        inverseDistortionCoefficients = None
        mappedDistortion              = False

    # Fetch parameters and convert
    
    pixelSize             = float(sim["CCD/PixelSize"])
    focalLength           = float(sim["Camera/FocalLength/ConstantValue"]) * 1000.0  # [m] -> [mm]
    raPlatform            = np.deg2rad(float(sim["ObservingParameters/RApointing"]))
    decPlatform           = np.deg2rad(float(sim["ObservingParameters/DecPointing"]))
    solarPanelOrientation = np.deg2rad(float(sim["Platform/SolarPanelOrientation"]))
    focalPlaneAngle       = np.deg2rad(float(sim["Camera/FocalPlaneOrientation/ConstantValue"]))

    telescopeGroup       = sim["Telescope/GroupID"]
    if telescopeGroup == "Custom":
        azimuthTelescope = np.deg2rad(float(sim["Telescope/AzimuthAngle"]))
        tiltTelescope    = np.deg2rad(float(sim["Telescope/TiltAngle"]))
    elif telescopeGroup == "Fast":
        azimuthTelescope = np.deg2rad(float(sim["CameraGroups/AzimuthAngle"][4]))
        tiltTelescope    = np.deg2rad(float(sim["CameraGroups/TiltAngle"][4]))
    else:
        idx = int(telescopeGroup)-1
        azimuthTelescope = np.deg2rad(float(sim["CameraGroups/AzimuthAngle"][idx]))
        tiltTelescope    = np.deg2rad(float(sim["CameraGroups/TiltAngle"][idx]))

    ccdZeroPointX = CCD[ccdCode]['zeroPointXmm']
    ccdZeroPointY = CCD[ccdCode]['zeroPointYmm']
    ccdAngle      = CCD[ccdCode]['angle']

    # Get the focal plane coordinates

    xFPmm, yFPmm = pixelToFocalPlaneCoordinates(xCCDpix, yCCDpix, pixelSize,
                                                ccdZeroPointX, ccdZeroPointY, ccdAngle)

    # If required, undistort them

    if includeFieldDistortion:
        if mappedDistortion:
            xFPmm, yFPmm = mappedDistortedToUndistortedFocalPlaneCoordinates(xFPmm, yFPmm,
                                                                             pathToPsfFile)
        else:
            xFPmm, yFPmm = distortedToUndistortedFocalPlaneCoordinates(xFPmm, yFPmm,
                                                                       inverseDistortionCoefficients,
                                                                       focalLength)

    # Get the corresponding sky coordinates

    ra, dec = focalPlaneToSkyCoordinates(xFPmm, yFPmm,
                                         raPlatform, decPlatform, solarPanelOrientation,
                                         tiltTelescope, azimuthTelescope,
                                         focalPlaneAngle, focalLength)

    # That's it!
    
    return ra, dec





def perturbPlatformPointing(x, y, z, ra, dec):

    """Perturbation to platform pointing

    This function determine a small perturbation in the platform's pointing
    given the euler angles (yaw, pitch, roll -> x, y, z) and the equatorial
    angles (ra, dec). I.e. a calculation from relative the euler angles to
    the relative equatorial coordinates.

    Parameters
    ----------
    x : float
        Euler angle perturbation in yaw [rad]
    y : float
        Euler angle perturbation in pitch [rad]
    z : float
        Euler angle perturbation in roll [rad]
    ra : float
        Right ascension of platform pointing [rad]
    dec : float
        Declination of platform pointing [rad]

    Return
    ------
    Perturbation of (RA, Dec).
    """

    # Rotation matrix for euler angles
    
    R = np.array([[ 0, -z,  y],
                  [ z,  0, -x],
                  [-y,  x,  0]])

    # Equatorial rotation matrix
    
    A = np.array([[np.cos(ra)*np.sin(dec)],
                  [np.sin(ra)*np.sin(dec)],
                  [1]])

    # That's it!
    
    return np.dot(R,A).T





# def matrixMisalignment(x, y, z):

#     """TODO not used yet
#     Parameters
#     ----------

#     Return
#     ------
#     """

#     r11 = + np.cos(x)*np.cos(z) - np.sin(x)*np.sin(z)*np.sin(y)
#     r12 = - np.cos(x)*np.sin(z) - np.sin(x)*np.cos(z)*np.cos(y)
#     r13 = + np.sin(x)*np.sin(z)
#     r21 = + np.sin(x)*np.cos(z) + np.cos(x)*np.sin(z)*np.cos(y)
#     r22 = - np.sin(x)*np.sin(z) - np.cos(x)*np.cos(z)*np.cos(y)
#     r23 = - np.cos(x)*np.sin(z)
#     r31 = + np.sin(z)*np.sin(y)
#     r32 = + np.cos(z)*np.sin(y)
#     r33 = - np.cos(y)
    
#     R = np.array([[r11, r12, r13],
#                   [r21, r22, r23],
#                   [r31, r32, r33]])

#     # 
#     return R
