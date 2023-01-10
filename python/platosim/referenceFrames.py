import h5py
import math
import os
import numpy as np

from numpy import *
from numpy.polynomial import Polynomial
from numpy.linalg import norm





# CCD configuration
#
# ccdCode:      1, 2, 3, 4: nominal cameras; 1F, 2F, 3F, 4F: fast cameras
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
    '1'  : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 0,    'zeroPointXmm':  -1.3, 'zeroPointYmm': +82.48, 'angle': pi},
    '2'  : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 0,    'zeroPointXmm':  -1.3, 'zeroPointYmm': +82.48, 'angle': 3*pi/2},
    '3'  : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 0,    'zeroPointXmm':  -1.3, 'zeroPointYmm': +82.48, 'angle': 0},
    '4'  : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 0,    'zeroPointXmm':  -1.3, 'zeroPointYmm': +82.48, 'angle': pi/2},
    '1F' : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 2255, 'zeroPointXmm':  -1.3, 'zeroPointYmm': +82.48, 'angle': pi},
    '2F' : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 2255, 'zeroPointXmm':  -1.3, 'zeroPointYmm': +82.48, 'angle': 3*pi/2},
    '3F' : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 2255, 'zeroPointXmm':  -1.3, 'zeroPointYmm': +82.48, 'angle': 0},
    '4F' : {'Nrows': 4510, 'Ncols': 4510, 'firstRow': 2255, 'zeroPointXmm':  -1.3, 'zeroPointYmm': +82.48, 'angle': pi/2}
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











def sunSkyCoordinatesAwayfromPlatformPointing(raPlatform, decPlatform, solarPanelOrientation=0.0):

    """
    Derive the location of the Sun which we assume to always be 180 degrees away from the platform pointing
    in the middle of the total time series.

    INPUT: raPlatform:  right ascension of the pointing of the Platform [rad]
           decPlatform: declination of the pointing of the Platform     [rad]
           solarPanelOrientation: orientation of the solar panel        [rad]
                                  (0, pi/2, pi, 3pi/2) for quarters (Q1,Q2,Q3,Q4)

    OUTPUT: raSun:  right ascension of the sun [rad]
            decSun: declination of the sun [rad]
    """

    lambdaPlatform, betaPlatform = equatorial2ecliptic(raPlatform, decPlatform)
    lambdaSun = lambdaPlatform - np.pi + solarPanelOrientation
    if (lambdaSun < 0.0): lambdaSun += 2.0 * np.pi
    raSun, decSun = ecliptic2equatorial(lambdaSun, 0.0)

    return raSun, decSun












def skyToPlatformCoordinates(raStar, decStar, raPlatform, decPlatform, solarPanelOrientation):

    """
    PURPOSE: Convert the equatorial sky coordinates (alpha, delta) of a star to cartesian platform coordinates
             (xSC, ySC, zSC)

    INPUT: raStar:                 right ascension of the star                               [rad]
           decStar:                declination of the star                                   [rad]
           raPlatform:             right ascension of the platform roll axis                 [rad]
           decPlatform:            declination of the platform roll axis                     [rad]
           solarPanelOrientation:  (0,pi/2,pi,3pi/2) for quarters (Q1,Q2,Q3,Q4)              [rad]

    OUTPUT: xSC, ySC, zSC: normalized cartesian coordinates of the direction of the star in the spacecraft reference frame.

    REMARK: Reference documents: PLATO-KUL-PL-TN-0001

    """

    # Get the sky position of the Sun (ra, dec) [rad]

    raSun, decSun = sunSkyCoordinatesAwayfromPlatformPointing(raPlatform, decPlatform, solarPanelOrientation)

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

    # Compute the cartesian coordinates of the star in the equatorial reference frame

    starEQ = np.array([cos(decStar)*cos(raStar), cos(decStar)*sin(raStar), sin(decStar)])

    # Transform these coordinates to the corresponding ones in the focal plane reference frame:

    starSC = np.dot(rotEQ2SC, starEQ)

    # That's it

    return starSC[0], starSC[1], starSC[2]









def platformToSkyCoordinates(xSC, ySC, zSC, raPlatform, decPlatform, solarPanelOrientation):

    """
    PURPOSE: Convert the cartesian platform coordinates (xSC, ySC, zSC) of a point to equatorial
             sky coordinates (alpha, delta). The units of the platform coordinates are arbitrary,
             since the output of this function is only a sky _direction_.

    INPUT: xSC:                    x-coordinate in the spacecraft reference frame  [arbitrary]
           ySC:                    y-coordinate in the spacecraft reference frame  [same unit as xSC]
           zSC:                    z-coordinate in the spacecraft reference frame  [same unit as xSC]
           raPlatform:             right ascension of the platform roll axis       [rad]
           decPlatform:            declination of the platform roll axis           [rad]
           solarPanelOrientation:  (0,pi/2,pi,3pi/2) for quarters (Q1,Q2,Q3,Q4)    [rad]

    OUTPUT: Equatorial coordinates of the direction of the vector (xSC, ySC, zSC)  [rad]

    REMARK: Reference documents: PLATO-KUL-PL-TN-0001

    """

    vecSC = np.array([xSC, ySC, zSC])

    # Get the sky position of the Sun (ra, dec) [rad]

    raSun, decSun = sunSkyCoordinatesAwayfromPlatformPointing(raPlatform, decPlatform, solarPanelOrientation)

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


    # Transform the unnormalized focal plane coordinates to the corresponding ones in the equatorial reference frame

    vecEQ = np.dot(rotSC2EQ, vecSC)

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












def skyToFocalPlaneCoordinates(raStar, decStar, raPlatform, decPlatform, solarPanelOrientation, tiltAngle, azimuthAngle, focalPlaneAngle, focalLength):

    """
    PURPOSE: Convert the equatorial sky coordinates (alpha, delta) of a star to undistorted normalized focal plane coordinates (xFP, yFP),
             assuming a spherical pinhole camera.

    INPUT: raStar:                 right ascension of the star                               [rad]
           decStar:                declination of the star                                   [rad]
           raPlatform:             right ascension of the platform roll axis                 [rad]
           decPlatform:            declination of the platform roll axis                     [rad]
           solarPanelOrientation:  (0,pi/2,pi,3pi/2) for quarters (Q1,Q2,Q3,Q4)              [rad]
           tiltAngle:              tilt angle of the telescope w.r.t. platform z-axis        [rad]
           azimuthAngle:           azimuth angle of the telescope on the platform            [rad]
           focalPlaneAngle:        angle between the Y_TL axis and the Y_FP axis: gamma_FP   [rad]
           focalLength:            focal length of the camera. Unit: see remarks.

    OUTPUT: xFP, yFP: normalized cartesian coordinates of the project star in the focal plane in the FP reference frame.
                      Unit: see remarks:

    REMARK: - The unit of the cartesian focal plane coordinates is the same as the one of focalLength.
              focalLength can be expressed in e.g. [mm] or in [pixels]. If focalLength == 1.0, then the corresponding
              focal plane coordinates are called "normalized coordinates."
            - Reference documents: PLATO-KUL-PL-TN-0001 and PLATO-DLR-PL-TN-016

    """

    # Get the sky position of the Sun (ra, dec) [rad]

    raSun, decSun = sunSkyCoordinatesAwayfromPlatformPointing(raPlatform, decPlatform, solarPanelOrientation)

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

    rotAzimuth = np.array([[ cos(azimuthAngle), sin(azimuthAngle), 0],   \
                           [-sin(azimuthAngle), cos(azimuthAngle), 0],   \
                           [        0        ,          0,         1]])

    rotTilt = np.array([[cos(tiltAngle), 0, -sin(tiltAngle)], \
                        [     0        , 1,        0       ], \
                        [sin(tiltAngle), 0,  cos(tiltAngle)]])

    rotSC2TL = np.dot(rotAzimuth.T, np.dot(rotTilt, rotAzimuth))

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

    xFPmm = -focalLength * starFP[0]/starFP[2]
    yFPmm = -focalLength * starFP[1]/starFP[2]

    # That's it

    return xFPmm, yFPmm












def focalPlaneToSkyCoordinates(xFP, yFP, raPlatform, decPlatform, solarPanelOrientation, tiltAngle, azimuthAngle, focalPlaneAngle, focalLength):

    """
    PURPOSE: Convert the undistorted normalized focal plane coordinates (xFP, yFP) of a star to equatorial sky coordinates (alpha, delta),
             assuming a spherical pinhole camera.

    INPUT: xFP:             undistorted normalized cartesian x-coordinate in the focal plane reference frame [same unit as focalLength]
           yFP:             undistorted normalized cartesian y-coordinate in the focal plane reference frame [same unit as focalLength]
           raPlatform:      right ascension of the platform pointing axis             [rad]
           decPlatform:     declination of the platform pointing axis                 [rad]
           solarPanelOrientation: (0,pi/2,pi,3pi/2) for quarters (Q1,Q2,Q3,Q4)        [rad]
           tiltAngle:       tilt angle of the telescope w.r.t. platform z-axis        [rad]
           azimuthAngle:    azimuth angle of the telescope on the platform            [rad]
           focalPlaneAngle: angle between the Y_TL axis and the Y_FP axis: gamma_FP   [rad]
           focalLength:     focal length of the camera.                               [mm]

    OUTPUT: raStar, decStar: Equatorial sky coordinates, right ascension and declination, of the star [rad]

    REMARK: The transformation assumes that the pinhole reverses the image.
    """

    # Get the sky position of the Sun (ra, dec) [rad]

    raSun, decSun = sunSkyCoordinatesAwayfromPlatformPointing(raPlatform, decPlatform, solarPanelOrientation)

    # Undo the reverse-image projection effect of the pinhole

    vecFP = np.array([-xFP/focalLength, -yFP/focalLength, 1.0], dtype=object)

    # Compute the rotation matrix to convert cartesian coordinates in the focal plane reference frame to
    # cartesian coordinates in the telescope reference frame

    rotFP2TL = np.array([[cos(focalPlaneAngle), -sin(focalPlaneAngle), 0],  \
                         [sin(focalPlaneAngle),  cos(focalPlaneAngle), 0],  \
                         [          0          ,           0         , 1]])

    # Compute the rotation matrix to convert cartesian coordinates in the telescope reference frame to
    # cartesian coordinates in the spacecraft reference frame

    rotAzimuth = np.array([[cos(azimuthAngle), -sin(azimuthAngle), 0],   \
                           [sin(azimuthAngle),  cos(azimuthAngle), 0],   \
                           [        0        ,          0,         1]])

    rotTilt = np.array([[ cos(tiltAngle), 0, sin(tiltAngle)], \
                        [     0         , 1,        0      ], \
                        [-sin(tiltAngle), 0, cos(tiltAngle)]])

    rotTL2SC = np.dot(rotAzimuth, np.dot(rotTilt, rotAzimuth.T))


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












def pixelCoordinates2FocalPlaneAngles(xCCD, yCCD, ccdCode, pixelSize, focalLength):

    """
    PURPOSE: Given the real-valued CCD pixel coordinates, compute the location angles of the star in the focal plane. These 
             calculations are based on the custom CCD position and orientation angle (so not from file!).

    INPUT:
        xCCD        : real-valued x-coordinate on the CCD (column)  [pix]
        yCCD        : real-valued y-coordinates on the CCD (row)    [pix]
        ccdCode     : one of: '1','2','3','4','1F','2F','3F','4F', depending on normal of fast cam
        pixelSize   : size of one square pixel   [microns]
        focalLength : focal length of the camera [mm]


    OUTPUT:
        angleFromOpticalAxis : angular distance from the optical axis               [rad]
        azimuthFromXAxis     : azimuth angle from the X-axis of the reference CCD   [rad]
    """

    ccdZeroPointX = CCD[ccdCode]["zeroPointXmm"]
    ccdZeroPointY = CCD[ccdCode]["zeroPointYmm"]
    ccdAngle      = CCD[ccdCode]["angle"]

    xmm, ymm = pixelToFocalPlaneCoordinates(xCCD, yCCD, pixelSize, ccdZeroPointX, ccdZeroPointY, ccdAngle)
    angleFromOpticalAxis = gnomonicRadialDistanceFromOpticalAxis(xmm,ymm,focalLength)
    azimuthFromXAxis = np.arctan2(ymm,xmm)

    return angleFromOpticalAxis, azimuthFromXAxis












def focalPlaneAngles2pixelCoordinates(angleFromOpticalAxis, azimuthFromXAxis, ccdCode, pixelSize, focalLength):

    """
    PURPOSE: given the location angles of the star in the focal plane, compute the real-valued CCD pixel coordinates. These 
             calculations are based on the custom CCD position and orientation angle (so not from file!).

    INPUT:
        angleFromOpticalAxis : angular distance from the optical axis             [rad]
        azimuthFromXAxis     : azimuth angle from the X-axis of the reference CCD [rad]
        ccdCode              : one of: '1','2','3','4','1F','2F','3F','4F', depending on normal of fast cam
        pixelSize            : size of one square pixel   [microns]
        focalLength          : focal length of the camera [mm]

    OUTPUT:
        xCCD : real-valued x-coordinate on the CCD (column)  [pix]
        yCCD : real-valued y-coordinates on the CCD (row)    [pix]
    """

    ccdZeroPointX = CCD[ccdCode]["zeroPointXmm"]
    ccdZeroPointY = CCD[ccdCode]["zeroPointYmm"]
    ccdAngle      = CCD[ccdCode]["angle"]

    xFPmm, yFPmm = focalPlaneCoordinatesFromGnomonicRadialDistance(angleFromOpticalAxis, focalLength, inPlaneRotation=azimuthFromXAxis)
    xCCD, yCCD = focalPlaneToPixelCoordinates(xFPmm, yFPmm, pixelSize, ccdZeroPointX, ccdZeroPointY, ccdAngle)

    return xCCD, yCCD









def telescopeToUndistortedFocalPlaneCoordinates(xTL, yTL, zTL, focalLength, focalPlaneAngle):

    """
    PURPOSE: given cartesian coordinates in the telescope reference frame, compute the undistorted
             coordinates in the focal plane reference frame.

    INPUT: xTL:             x-coordinate in the focal plane reference frame      [arbitrary unit]
           yTL:             y-coordinate in the focal plane reference frame      [  "         " ]
           zTL:             z-coordinate in the focal plane reference frame      [  "         " ]
           focalLength:     Focal length of the camera.                          [mm]
           focalPlaneAngle: focal plane orientation angle                        [rad]

    OUTPUT: xFP, yFP: cartesian coordinates in the focal plane reference frame.  [mm]
    """

    rotTL2FP = np.array([[ cos(focalPlaneAngle), sin(focalPlaneAngle), 0],  \
                         [-sin(focalPlaneAngle), cos(focalPlaneAngle), 0],  \
                         [          0          ,           0         , 1]])

    vecTL = np.array([xTL, yTL, zTL])
    vecFP = np.dot(rotTL2FP, vecTL)

    # Convert the units to the one of focalLength (usually [mm]), and normalize the coordinates
    # to take into account the pinhole camera projection.

    xFPmm = -focalLength * vecFP[0]/vecFP[2]
    yFPmm = -focalLength * vecFP[1]/vecFP[2]

    return xFPmm, yFPmm










def undistortedFocalPlaneToTelescopeCoordinates(xFP, yFP, focalLength, focalPlaneAngle):

    """
    PURPOSE: given cartesian coordinates in the focal plane reference frame, compute the
             coordinates in the telescope reference frame.

    INPUT: xFP:             x-coordinate in the focal plane reference frame   [mm]
           yFP:             y-coordinate in the focal plane reference frame   [mm]
           focalLength:     Focal length of the camera.                       [mm]
           focalPlaneAngle: focal plane orientation angle                     [rad]

    OUTPUT: xTL, yTL, zTL: cartesian coordinates in the telescope reference frame.

    REMARK: Because of the projection degeneracy, the telescope coordinates are computed on the unit sphere
    """


    vecFP = np.array([-xFP/focalLength, -yFP/focalLength, 1.0])

    # Compute the rotation matrix to convert cartesian coordinates in the focal plane reference frame to
    # cartesian coordinates in the telescope reference frame

    rotFP2TL = np.array([[cos(focalPlaneAngle), -sin(focalPlaneAngle), 0],  \
                         [sin(focalPlaneAngle),  cos(focalPlaneAngle), 0],  \
                         [          0          ,           0         , 1]])

    vecTL = np.dot(rotFP2TL, vecFP)
    vecTL /= norm(vecTL)

    return vecTL[0], vecTL[1], vecTL[2]












def undistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm, distortionCoefficients, focalLength):

    """
    PURPOSE:      Convert from undistorted to distorted normalized focal plane coordinates using the analytic distortion model.

    INPUTS:       xFPmm  undistorted normalized focal plane x-coordinate [mm]
                  yFPmm  undistorted normalized focal plane y-coordinate [mm]
                  distortionCoefficients  List of polynomial coefficients
                  focalLength: Focal length [mm]

    OUTPUTS:      (xFPdist, yFPdist) distorted x and y coordinates [mm]

    REMARK:       This is not the prefered method for detectors with mapped PSF, since these use a mapped distortion model.
                  For such models use the function mappedUndistortedToDistortedFocalPlaneCoordinates.

    Note: Example of distortion coefficients: [0.32419,  0.0232909,  0.407979, 0.00022463, 0.000217599, 0.000381958, 0.000963902]
    """

    radialCoefficients = [0, 0, 0, distortionCoefficients[0], 0, distortionCoefficients[1], 0, distortionCoefficients[2]]
    distortionPolynomial = Polynomial(radialCoefficients)

    angle = arctan2(yFPmm, xFPmm)    # Position angle on the focal plane [radians]

    rFP = sqrt(xFPmm**2 + yFPmm**2) / focalLength              # Undistorted radial distance [normalised pixels]
    radialDistortion = distortionPolynomial(rFP) * focalLength
    tangentialDistortion = rFP**2 * (distortionCoefficients[5] * cos(angle) + distortionCoefficients[6] * sin(angle)) * focalLength

    xFPdist = xFPmm + cos(angle) * (radialDistortion + tangentialDistortion) + distortionCoefficients[3] * rFP**2 * focalLength
    yFPdist = yFPmm + sin(angle) * (radialDistortion + tangentialDistortion) + distortionCoefficients[4] * rFP**2 * focalLength

    return xFPdist, yFPdist











def mappedUndistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm, pathToPsfFile):

    """
    PURPOSE:      Convert from undistorted to distorted normalized focal plane coordinates using a mapped distortion model.

    INPUTS:       xFPmm:  undistorted normalized focal plane x-coordinate [mm]
                  yFPmm:  undistorted normalized focal plane y-coordinate [mm]w
                  pathToFoPsfFile: absolute path to the PSF file that the Coordinates map

    OUTPUTS:      (xFPdist, yFPdist) distorted x and y coordinates [mm]

    REMARK:       This is the prefered method for detectors with mapped PSF.
                  For detectors that use analytic PSF use the function undistortedToDistortedFocalPlaneCoordinates.

    """

    # Check the path to the PSF file excists

    if not os.path.exists(pathToPsfFile):
        if os.path.exists(os.environ["PLATO_PROJECT_HOME"] + "/" + pathToPsfFile):
            pathToPsfFile = os.environ["PLATO_PROJECT_HOME"] + "/" + pathToPsfFile
        else:
            print("Error: {} is not a valid path name for mapped PSF".format(pathToPsfFile))

    # We open the psf file where the coordinate transformation matrix should be in. If the matrix isn't in the file, we raise an error.

    psfFile = h5py.File(pathToPsfFile, "r")
    if not "Coordinates map" in psfFile.keys():
        print("Error: No transformation map given in psf file, mapped distortion is not possible.")
        return
    else:
        coordMap = psfFile["Coordinates map"]


    # Calculate the distance of the undistorted coordinates with the input coordinates.

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

    distanceFromPoint = [max( (xFPmm - x)**2, (yFPmm - y)**2) for x, y in zip(xUndis, yUndis)]

    # We should select the closest four undistorted points to the input point

    idx = list(range(len(distanceFromPoint)))
    idx.sort(key=lambda i: distanceFromPoint[i])
    idx = [ idx[i] for i in [0, 1, 2, 3]]


    # We sort these points
    first_points_idx  = [i for i in idx if (yUndis[i] < yFPmm)]
    first_points_idx.sort(key=lambda idx: xUndis[idx])
    second_points_idx = [i for i in idx if (yUndis[i] >= yFPmm)]
    second_points_idx.sort(key=lambda idx: xUndis[idx])

    idx = first_points_idx + second_points_idx


    closestX, closestY = np.array([xUndis[idx[0]], xUndis[idx[1]], xUndis[idx[2]], xUndis[idx[3]]]), np.array([yUndis[idx[0]], yUndis[idx[1]], yUndis[idx[2]], yUndis[idx[3]]])

    # We can write the points (xFPmm, yFPmm) as a linear combination of the
    # four closests points around this point.

    oPointIdx = [3, 2, 1, 0]
    area = (closestX[0] - closestX[3])*(closestY[0]-closestY[3])

    constants = [abs( (closestX[oPointIdx[i]] - xFPmm) *
                      (closestY[oPointIdx[i]] - yFPmm))/ area
                 for i in np.arange(4)]



    closestXdist = np.array([ xDist[idx[i]] for i in [0,1,2,3]])
    closestYdist = np.array([ yDist[idx[i]] for i in [0,1,2,3]])


    xFPdist = sum([constants[i]*closestXdist[i] for i in np.arange(4)])
    yFPdist = sum([constants[i]*closestYdist[i] for i in np.arange(4)])

    return round(xFPdist, 4), round(yFPdist, 4)










def distortedToUndistortedFocalPlaneCoordinates(xFPdist, yFPdist, inverseDistortionCoefficients, focalLength):

    """
    PURPOSE:     Convert from distorted to undistorted normalized focal plane coordinates using the analytic distortion model.

    INPUTS:      xFPdist  Distorted normalized focal plane x-coordinate [mm]
                 yFPdist  DIstorted normalized focal plane y-coordinate [mm]
                 inverseDistortionCoefficients  List of polynomial coefficients
                 focalLength: Focal length [mm]

    OUTPUTS:     (xFPmm, yFPmm) distorted x and y coordinates [mm]

    REMARK:     This is not the prefered method for detectors with mapped PSF, since these use a mapped distortion model.
                For such models use the function mappedDistortedToUndistortedFocalPlaneCoordinates.

    Note: Example of inverse distortion coefficients: [-0.323487, 0.268344, -0.435473, -0.00019304, -0.000176961, -0.000321713, -0.000827654]
    """

    inverseCoefficientsRadial = [0, 0, 0, inverseDistortionCoefficients[0], 0, inverseDistortionCoefficients[1], 0, inverseDistortionCoefficients[2]]
    inverseDistortionPolynomialRadial = Polynomial(inverseCoefficientsRadial)

    angle = arctan2(yFPdist, xFPdist)     # Position angle on the focal plane [radians]

    rFP = sqrt(xFPdist**2 + yFPdist**2) / focalLength                   # Distorted radial distance [normalised pixels]
    radialDistortion = inverseDistortionPolynomialRadial(rFP) * focalLength         # Distortion [mm] -> negative!
    tangentialDistortion = rFP**2 * (inverseDistortionCoefficients[5] * cos(angle) + inverseDistortionCoefficients[6] * sin(angle)) * focalLength

    xFPmm = xFPdist + cos(angle) * (radialDistortion + tangentialDistortion) + inverseDistortionCoefficients[3] * rFP**2 * focalLength
    yFPmm = yFPdist + sin(angle) * (radialDistortion + tangentialDistortion) + inverseDistortionCoefficients[4] * rFP**2 * focalLength

    return xFPmm, yFPmm











def mappedDistortedToUndistortedFocalPlaneCoordinates(xFPdist, yFPdist, pathToPsfFile):

    """
    PURPOSE:     Convert from distorted to undistorted normalized focal plane coordinates using the mapped distortion model.

    INPUTS:      xFPdist  Distorted normalized focal plane x-coordinate [mm]
                 yFPdist  DIstorted normalized focal plane y-coordinate [mm]
                 pathToFoPsfFile: absolute path to the PSF file that the Coordinates map

    OUTPUTS:     (xFPmm, yFPmm) distorted x and y coordinates [mm]

    REMARK:     This is the prefered method for detectors with mapped PSF, since these use a mapped distortion model.
                For detectors that use analytic PSF use the function distortedToUndistortedFocalPlaneCoordinates.
    """

    # Check the path to the PSF file excists
    if not os.path.exists(pathToPsfFile):
        if os.path.exists(os.environ["PLATO_PROJECT_HOME"] + "/" + pathToPsfFile):
            pathToPsfFile = os.environ["PLATO_PROJECT_HOME"] + "/" + pathToPsfFile
        else:
            print("Error: {} is not a valid path name for mapped PSF".format(pathToPsfFile))

    # We open the psf file where the coordinate transformation matrix should be in. If this map isn't in the file, we raise an error.
    psfFile = h5py.File(pathToPsfFile, "r")
    if not "Coordinates map" in psfFile.keys():
        print("Error: No transformation map given in psf file, mapped distortion is not possible.")
        return
    else:
        coordMap = psfFile["Coordinates map"]

    # Calculate the distance of the distorted coordinates with the distorted input coordinates.

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

    distanceFromPoint = [max( (xFPdist - x)**2, (yFPdist - y)**2) for x, y in zip(xDist, yDist)]

    # We should select the four closest distorted points to the input point

    idx = list(range(len(distanceFromPoint)))
    idx.sort(key=lambda i: distanceFromPoint[i])
    idx = [idx[i] for i in [0, 1, 2, 3]]

    # We sort these points
    first_points_idx = [i for i in idx if (yDist[i] < yFPdist)]
    first_points_idx.sort(key=lambda idx: xDist[idx])
    second_points_idx = [i for i in idx if (yDist[i] >= yFPdist)]
    second_points_idx.sort(key=lambda idx: xDist[idx])

    idx = first_points_idx + second_points_idx

    closestX, closestY = np.array([xDist[idx[0]], xDist[idx[1]], xDist[idx[2]], xDist[idx[3]]]), np.array([yDist[idx[0]], yDist[idx[1]], yDist[idx[2]], yDist[idx[3]]])

    # We can write the points (xFPmm, yFPmm) as a linear combination of the
    # four closest points around this point.

    oPointIdx = [3, 2, 1, 0]
    area = (closestX[0] - closestX[3])*(closestY[0] - closestY[3])

    constants = [abs( (closestX[oPointIdx[i]] - xFPdist) *
                      (closestY[oPointIdx[i]] - yFPdist))
                 for i in np.arange(4)]
    constants = [ constant / sum(constants) for constant in constants]

    closestXund = np.array([ xUndis[idx[i]] for i in [0,1,2,3]])
    closestYund = np.array([ yUndis[idx[i]] for i in [0,1,2,3]])


    xFPmm = sum([constants[i]*closestXund[i] for i in np.arange(4)])
    yFPmm = sum([constants[i]*closestYund[i] for i in np.arange(4)])

    return round(xFPmm, 4), round(yFPmm, 4)











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
           focalLength focal length of the camera [mm]

    OUTPUT: the angular distance of the star w.r.t. the optical axis [rad]

    """

    tanx = xFP / focalLength
    tany = yFP / focalLength

    angularDistance = np.arccos(1.0/sqrt(1.0 + tanx*tanx + tany*tany));

    # Take care that the angle is between [0, 2*PI]

    # could be simplified into "return angularDistance % (2*np.pi)"
    if angularDistance < 0.0:
        angularDistance += 2.0 * np.pi
    elif angularDistance > 2.0 * np.pi:
        angularDistance -= 2.0 * np.pi

    # That's it!

    return angularDistance;









def focalPlaneCoordinatesFromGnomonicRadialDistance(angularDistance, focalLength, inPlaneRotation=0.):

    """
    Calculate the xFP,yFP focal plane coordinates from the gnomonic
    radial distance with respect to the optical axis in the focal plane

    INPUT: angularDistance: angular distance of the star w.r.t. the optical axis [rad]
           focalLength    : focal length of the camera                           [mm]
           inPlaneRotation: angle from the xFP axis to the target (default=0)    [rad]

    OUTPUT: xFP  Focal plane x-coordinate [mm]
            yFP  Focal plane y-coordinate [mm]

    """

    D = focalLength * np.tan(angularDistance)

    xFP = D * np.cos(inPlaneRotation)
    yFP = D * np.sin(inPlaneRotation)

    return xFP,yFP












def computeCCDcornersInFocalPlane(ccdCode, pixelSize):

    """
    PURPOSE: Get the (x,y) coordinates of each of the 4 corners of the exposed part of the CCD
             in the FP' reference system.  These calculations are based on the custom 
             CCD position and orientation angle (so not from file!).

    INPUT: ccdCode:   one of the following: '1', '2', '3', '4', '1F', '2F', '3F', '4F'
           pixelSize: size of 1 pixel (micron)

    OUTPUT: cornersXmm: x-coordinates of each of the corners in the FP' reference system [mm]
            cornersYmm: y-coordinates of each of the corners in the FP' reference system [mm]
    """

    # Get the pixel coordinates of the 4 corners of the exposed part of the CCD
    # Note that the x-direction corresponds to the CCD columns, and the y-direction to the CCD rows.

    Nrows = CCD[ccdCode]["Nrows"]
    Ncols = CCD[ccdCode]["Ncols"]
    firstRow = CCD[ccdCode]["firstRow"]

    cornersXpix = array([0.0, Ncols, Ncols, 0.0])
    cornersYpix = array([firstRow, firstRow, Nrows, Nrows])

    # Convert to the x,y coordinates in the FP' reference frame

    zeroPointXmm = CCD[ccdCode]["zeroPointXmm"]
    zeroPointYmm = CCD[ccdCode]["zeroPointYmm"]
    ccdAngle     = CCD[ccdCode]["angle"]

    cornersXmm, cornersYmm = pixelToFocalPlaneCoordinates(cornersXpix, cornersYpix, pixelSize, zeroPointXmm, zeroPointYmm, ccdAngle)

    # That's it

    return cornersXmm, cornersYmm







def getCCDandPixelCoordinates(raStar, decStar, raPlatform, decPlatform, solarPanelOrientation,
                              tiltAngle, azimuthAngle, focalPlaneAngle, focalLength, pixelSize,
                              includeFieldDistortion, normal, mappedDistortion=False,
                              distortionCoefficients=None, pathToPsfFile=None):

    """
    PURPOSE: Given the equatorial coordinates of a star, find out on which CCD it falls ('1', '2', ...)
             and compute the pixel coordinates of the star on this CCD. If the star doesn't fall on any of the CCDs
             then (None, None, None) is given as output.  These calculations use the custom CCD positions and 
             orientation angle (so not from file!).

    INPUT: raStar:                 right ascension of the star                               [rad]
           decStar:                declination of the star                                   [rad]
           raPlatform:             right ascension of the platform pointing axis             [rad]
           decPlatform:            declination of the platform pointing axis                 [rad]
           solarPanelOrientation:  (0,pi2/,pi,3pi/2) for quarters (Q1,Q2,Q3,Q4)              [rad]
           tiltAngle:              tilt angle of the telescope w.r.t. platform z-axis        [rad]
           azimuthAngle:           azimuth angle of the telescope on the platform            [rad]
           focalPlaneAngle:        angle between the Y_TL axis and the Y_FP axis: gamma_FP   [rad]
           focalLength:            focal length of the camera.                               [mm]
           pixelSize:              pixel size                                                [micron]
           includeFieldDistortion: True to include field distortion in coordinate transformations, false otherwise
           normal:                 True for the normal camera configuration, False for the fast cameras
           mappedDistortion:       True if we want mapped distortion (mapped from file psf) False if we have analytic psfs      
           distortionCoefficients: Coefficients of the polynomial describing the distortion for anlytic psf 
           pathToPsfFile         : Path to the PSF file for mapped PSFs to calculate mapped distortion


    OUTPUT: ccdCode: for normal camera: either '1', '2', '3', or '4'
                     for fast camera: either '1F', '2F', '3F', '4F'
                     if on no CCD: None
            xCCDpix: x-coordinate (column number) of the star on the CCD  [pix]
                     if on no CCD: None
            yCCDpix: y-coordinate (row number) of the star on the CCD  [pix]
                     if on no CCD: None

    Note: Example of distortion coefficients: [-0.0036696919678, 1.0008542317, -4.12553764817e-05, 5.7201219949e-06]
    """

    # Make sure that for the respective field distortion the proper information is given.
    
    if (includeFieldDistortion or includeFieldDistortion == "yes"):
        if (mappedDistortion and pathToPsfFile is None):
            print("Error: If mapped field distortion should be taken into account, a path to the psf file should be given")
            return
        elif ( (not mappedDistortion) and distortionCoefficients is None):
            print("Error: If analytic field distortion should be taken into account, the distortionCoefficients should be given")
            return
        
    # Select the proper CCD codes depending on whether we're dealing with the nominal or the fast cams

    if normal == True:
        ccdCodes = ['1', '2', '3', '4']
    else:
        ccdCodes = ['1F', '2F', '3F', '4F']


    # Compute the (x,y) coordinates in the FP reference system [mm]

    xFPmm, yFPmm = skyToFocalPlaneCoordinates(raStar, decStar, raPlatform, decPlatform, solarPanelOrientation,
                                              tiltAngle, azimuthAngle, focalPlaneAngle, focalLength)

    if (includeFieldDistortion == True) or (includeFieldDistortion == "yes"):
        if mappedDistortion:
            xFPmm, yFPmm = mappedUndistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm, pathToPsfFile)
        else: 
            xFPmm, yFPmm = undistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm, distortionCoefficients, focalLength)

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

    OUTPUT: raTelescope: right ascension of the optical axis of the telescope    [rad]
            decTelescope: declination of the optical axis of the telescope       [rad]

    """


    # Compute the rotation matrix to convert cartesian coordinates in the telescope reference frame to
    # cartesian coordinates in the spacecraft reference frame

    rotAzimuth = np.array([[cos(azimuthAngle), -sin(azimuthAngle), 0],   \
                           [sin(azimuthAngle),  cos(azimuthAngle), 0],   \
                           [        0        ,          0,         1]])

    rotTilt = np.array([[ cos(tiltAngle), 0, sin(tiltAngle)], \
                        [     0        , 1,        0       ], \
                        [-sin(tiltAngle), 0, cos(tiltAngle)]])

    rotTL2SC = np.dot(rotAzimuth, np.dot(rotTilt, rotAzimuth.T))


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











def getCameraGroupCoordinates(raPlatform, decPlatform, solarPanelOrientation=0):

    """
    PURPOSE: Calculate the ICRS coordinates of each camera group.

    INPUT: raPlatform:            Right Ascension of the platform pointing [deg]
           decPlatform:           Declination of the platform pointing [deg]
           solarPanelOrientation: Orientation of the solar panel [deg]
           
    OUTPUT: RA and Dec coordinates for each camera group.
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











def calculateSubfieldAroundCoordinates(subfieldSizeX, subfieldSizeY, raStar, decStar, raPlatform, decPlatform,
                                       solarPanelOrientation, tiltTelescope, azimuthTelescope,
                                       focalPlaneAngle, focalLength, pixelSize,
                                       includeFieldDistortion, normal,
                                       mappedDistortion=False, distortionCoefficients=None, pathToPsfFile=None ):

    """
    PURPOSE: Calculates the location of the subfield such that the star with coordinates (raStar, decStar)
             is centered in the subfield. The function also checks if there is enough space around the pixel
             feature/feature-447-CCD

    NOTE:    This function is used by setSubfieldAroundCoordinates() and usually does not need to be called by the user.

    INPUT: subfieldSizeX:          full width (# of columns) of the subfield                [pix]
           subfieldSizeY:          full height (#of rows) of the subfield                   [pix]
           raStar:                 right ascension of the star                              [rad]
           decStar:                declination of the star                                  [rad]
           raPlatform:             right ascension of the platform pointing axis            [rad]
           decPlatform:            declination of the platform pointing axis                [rad]
           solarPanelOrientation:  (0,pi/2,pi,3pi/2) for quarters (Q1,Q2,Q3,Q4)             [rad]
           tiltTelescope           tilt angle of the telescope w.r.t. platform z-axis       [rad]
           azimuthTelescope:       azimuth angle of the telescope on the platform           [rad]
           focalPlaneAngle:        angle between the Y_TL axis and the Y_FP axis: gamma_FP  [rad]
           focalLength:            focal length of the camera.                              [mm]
           pixelSize:              pixel size                                               [micron]
           includeFieldDistortion: True to include field distortion in coordinate transformations, false otherwise
           normal:                 True for the normal camera configuration, False for the fast cameras
           mappedDistortion:       True if we want mapped distortion (mapped from file psf) False if we have analytic psfs      
           distortionCoefficients: Coefficients of the polynomial describing the distortion for anlytic psf 
           pathToPsfFile:          Path to the PSF file for mapped PSFs to calculate mapped distortion       



    OUTPUT: ccdCode: "1", "2", "3" or "4" if nominal=True, "1F", "2F", "3F" or "4F" otherwise
            xCCDpix: x-coordinate of the star in pixels (i.e. column number)
            yCCDpix: y-coordinate of the star in pixels (i.e. row number)

    REMARKS: - If the coordinates do not fall on any CCD, an error message is shown
             - If the star is too close to the edge for the given subfield size, and error message is shown,
               followed by an exit(1)
             - Example of distortion coefficients: [-0.0036696919678, 1.0008542317, -4.12553764817e-05, 5.7201219949e-06]
    """

    # Find out that we have been given the correct distortion input parameters. If this is not the case raise error and return. 
    if (includeFieldDistortion or includeFieldDistortion == "yes"):
        if (mappedDistortion and pathToPsfFile is None):
            print("Error: If mapped field distortion should be taken into account, a path to the psf file should be given")
            return
        elif ( (not mappedDistortion) and distortionCoefficients is None):
            print("Error: If analytic field distortion should be taken into account, the distortionCoefficients should be given")
            return

    # Find out on which CCD the star falls, and the corresponding pixel coordinates

    ccdCode, xCCDpix, yCCDpix = getCCDandPixelCoordinates(raStar, decStar, raPlatform, decPlatform, solarPanelOrientation,
                                                          tiltTelescope, azimuthTelescope, focalPlaneAngle, focalLength,
                                                          pixelSize, includeFieldDistortion, normal,
                                                          mappedDistortion, distortionCoefficients, pathToPsfFile)

    # If the CCD code is None, the star does not fall on any ccd -> error

    if ccdCode == None:
        return None, None, None

    # If the star does fall on a CCD, check if it's not too close to the edge for the subfield to
    # be completely on the CCD.

    xCCDpix = int(xCCDpix)                # integer values
    yCCDpix = int(yCCDpix)
    firstRow = CCD[ccdCode]["firstRow"]     # different from nominal than for fast cams
    Ncols = CCD[ccdCode]["Ncols"]
    Nrows = CCD[ccdCode]["Nrows"]
    
    if     (xCCDpix - subfieldSizeX/2 < 0)        or (xCCDpix + subfieldSizeX/2 - 1 > Ncols-1)   \
        or (yCCDpix - subfieldSizeY/2 < firstRow) or (yCCDpix + subfieldSizeY/2 - 1 > Nrows-1):
        return None, None, None

    # That's it!

    return ccdCode, xCCDpix, yCCDpix












def skyToPixelCoordinates(sim, raStar, decStar, normal):
    """
    PURPOSE: Convert sky coordinates to pixel coordinates. These calculations are based on the custom 
             CCD position and orientation angle (so not from file!).
             
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

    if (sim["PSF/Model"] == "MappedFromFile"):
        includeFieldDistortion = True
        distortionCoefficients = None 
        pathToPsfFile          = sim["PSF/MappedFromFile/Filename"]
        mappedDistortion       = True
    elif (sim["Camera/IncludeFieldDistortion"] == "yes")  or (sim["Camera/IncludeFieldDistortion"] == "1") or (sim["Camera/IncludeFieldDistortion"]):
        distortionCoefficients = sim["Camera/FieldDistortion/ConstantCoefficients"]
        pathToPsfFile          = None
        mappedDistortion       = False
        includeFieldDistortion = True
    else:
        includeFieldDistortion = False
        distortionCoefficients = None
        pathToPsfFile          = None
        mappedDistortion       = False

        
    pixelSize             = float(sim["CCD/PixelSize"])
    focalLength           = float(sim["Camera/FocalLength/ConstantValue"]) * 1000.0                   # [m] -> [mm]
    raPlatform            = np.deg2rad(float(sim["ObservingParameters/RApointing"]))
    decPlatform           = np.deg2rad(float(sim["ObservingParameters/DecPointing"]))
    solarPanelOrientation = np.deg2rad(float(sim["Platform/SolarPanelOrientation"]))
    focalPlaneAngle       = np.deg2rad(float(sim["Camera/FocalPlaneOrientation/ConstantValue"]))
    azimuthTelescope      = np.deg2rad(float(sim["Telescope/AzimuthAngle"]))
    tiltTelescope         = np.deg2rad(float(sim["Telescope/TiltAngle"]))

    # Get the pixel coordinates on the CCD

    ccdCode, xCCDpixel, yCCDpixel = getCCDandPixelCoordinates(raStar, decStar, raPlatform, decPlatform, solarPanelOrientation,\
                                                              tiltTelescope, azimuthTelescope, focalPlaneAngle, focalLength, \
                                                              pixelSize, includeFieldDistortion, normal, mappedDistortion, \
                                                              distortionCoefficients, pathToPsfFile)

    return ccdCode, xCCDpixel, yCCDpixel















def pixelToSkyCoordinates(sim, ccdCode, xCCDpixel, yCCDpixel):

    """
    PURPOSE: Convert pixel coordinates to equatorial sky coordinates (ra,dec).  These calculations are based on the custom 
             CCD position and orientation angle (so not from file!).
             
    NOTE:   It is assumed that the configuration parameters in the sim object contains
            a correct (ra, dec) of the platform, a correct (azimuth, tilt) of the telescope,
            a valid value for the focal length, the pixel size, and that
            the switch to include distortion or not is set correctly.

    INPUT:  sim:        simulation for which the configuration file is adapted
            
            :    for nominal camera: either '1', '2', '3', '4'
                        for fast camera: either '1F', '2F', '3F', '4F'
            xCCDpixel:  x-coordinate (column-number) of the star on the CCD  [pix]
            yCCDpixel:  y-coordinate (row-number) of the star on the CCD     [pix]

    OUTPUT: raStar, decStar: Equatorial coordinates (right ascension and declination) of the star [rad]
    """

    if (sim["PSF/Model"] == "MappedFromFile"):
        includeFieldDistortion = True        
        inverseDistortionCoefficients = None 
        pathToPsfFile          = sim["PSF/MappedFromFile/Filename"]
        mappedDistortion       = True

    elif (sim["Camera/IncludeFieldDistortion"] == "yes")  or (sim["Camera/IncludeFieldDistortion"] == "1") or (sim["Camera/IncludeFieldDistortion"] == "True"):
            inverseDistortionCoefficients = sim["Camera/FieldDistortion/ConstantInverseCoefficients"]
            pathToPsfFile          = None
            mappedDistortion       = False
            includeFieldDistortion = True
           
    else:
        includeFieldDistortion        = False
        pathToPsfFile                 = None
        inverseDistortionCoefficients = None
        mappedDistortion              = False


    pixelSize             = float(sim["CCD/PixelSize"])
    focalLength           = float(sim["Camera/FocalLength/ConstantValue"]) * 1000.0                     # [m] -> [mm]
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

    xFPmm, yFPmm = pixelToFocalPlaneCoordinates(xCCDpixel, yCCDpixel, pixelSize, ccdZeroPointX, ccdZeroPointY, ccdAngle)

    # If required, undistort them

    if includeFieldDistortion:
        if mappedDistortion:
            xFPmm, yFPmm = mappedDistortedToUndistortedFocalPlaneCoordinates(xFPmm, yFPmm, pathToPsfFile)
        else:
            xFPmm, yFPmm = distortedToUndistortedFocalPlaneCoordinates(xFPmm, yFPmm, inverseDistortionCoefficients, focalLength)

    # Get the corresponding sky coordinates

    ra, dec = focalPlaneToSkyCoordinates(xFPmm, yFPmm, raPlatform, decPlatform, solarPanelOrientation, tiltTelescope, azimuthTelescope, focalPlaneAngle, focalLength)

    return ra, dec











def getPointingRepeatabilityError(ra, dec, rot, sigma=[3., 6.], quarter=[1, 8], show_table=False):
    """
    TODO under development!

    PURPOSE: 
             
    INPUT:

    OUTPUT:
    """

    # Coordinates
    ICRS = np.array([ra, dec, rot])

    # Pointing Reproducibility Error (PRE) in P/L reference frame (yaw, pitch, roll)
    # Here t stands for transverse direction and  
    sigma_t = sigma[0]/3600
    sigma_b = sigma[1]/3600

    # Find distribution within 3 sigma of req.
    tt = np.array([np.random.normal(0, t/sigma) for i in range(len(quarters))])
    bb = np.array([np.random.normal(0, b/sigma) for i in range(len(quarters))])

    # Corresponding yaw, pitch, roll
    y = tt
    z = 3 * y
    x = bb - z

    # ICRS pointing angles
    phi   = np.deg2rad(ra)
    theta = np.deg2rad(dec)

    # Find change to pointing for quarters
    coor = np.zeros((len(quarters), 4))
    for i in range(len(quarters)):
        data = changeOfPointing(x[i], y[i], z[i], phi, theta)[0]
        coor[i,:] = np.append(quarters[i], data)

    # Save file with relative pointing errors [deg]
    np.savetxt(f'{outdir}/PRE.txt', coor, fmt=['%i', '%0.8f', '%0.8f', '%0.8f'])



















