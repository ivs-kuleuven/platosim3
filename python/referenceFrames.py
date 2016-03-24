from numpy import *

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






def skyToAngularFocalPlaneCoordinates(raStar, decStar, raOpticalAxis, decOpticalAxis, focalPlaneAngle, plateScale, pixelSize):

    """
    PURPOSE: computes the (x,y) coordinates in the focal plane of a star with given equatorial coordinates

    INPUT: raStar:          right ascension of the star [rad]
           decStar:         declination of the star [rad]
           raOpticalAxis:   right ascension of the optical axis [rad]
           decOpticalAxis:  declination of the optical axis [rad]
           focalPlaneAngle: angle between the Y_FP axis and the Y'_FP axis: gamme_FP  [rad]
           plateScale:      [arcsec/micron]
           pixelSize:       [micrometer]

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









def angularToPlanarFocalPlaneCoordinates(xFPrad, yFPrad, focalLength):
    """
    PURPOSE: Convert from angular to planar focal plane coordinates, assuming no optical distortion.

    INPUT:   xFPrad   Angular focal plane x-coordinate [rad]
             yFPrad   Angular focal plane y-coordinate [rad]
    
    OUTPUT:  (xFPmm, yFPmm)    Planar focal plane x and y coordinates [mm]

    """

    xFPmm = tan(xFPrad) * focalLength
    yFPmm = tan(yFPrad) * focalLength

    return xFPmm, yFPmm







def inverseGnomonicProjectionFocalPlaneToSky(xFPprimeStar, yFPprimeStar, raOpticalAxis, decOpticalAxis, focalPlaneAngle, plateScale):

    """
    PURPOSE: computes the (x,y) coordinates in the focal plane of a star with given equatorial coordinates

    INPUT: xFPprimeStar:    x-coordinate of the projected star in the focal plane in the FP-prime system [mm]
           xFPprimeStar:    y-coordinate of the projected star in the focal plane in the FP-prime system [mm]
           raOpticalAxis:   right ascension of the optical axis [rad]
           decOpticalAxis:  declination of the optical axis [rad]
           focalPlaneAngle: angle between the Y_FP axis and the Y'_FP axis: gamme_FP  [rad]
           plateScale:      [arcsec/micron]

    OUTPUT: raStar:  right ascension of the star [rad]
            decStar: declination of the star [rad]
            
    REMARK: This function does not work when the star is located at 0, 0 (Zenit).
    """

    
    if isscalar(xFPprimeStar) and isscalar(yFPprimeStar) and xFPprimeStar == 0.0 and yFPprimeStar == 0.0:
        return raOpticalAxis, decOpticalAxis
    
    # Compute the conversion factor from [mm/radian] to [arcsec/pixel] and convert the
    # focal plane coordinates from [mm] to [rad].

    conversionFactor = 3600. * 180 / 1000 / pi / plateScale
    xFPprime = xFPprimeStar / conversionFactor
    yFPprime = yFPprimeStar / conversionFactor

    # Convert the FP' coordinates into FP coordinates

    xFP =  xFPprime * cos(focalPlaneAngle) - yFPprime * sin(focalPlaneAngle)
    yFP =  xFPprime * sin(focalPlaneAngle) + yFPprime * cos(focalPlaneAngle)

    # Project the focal plane in the "FP" coordinate system to the sky

    rho = sqrt(xFP*xFP+yFP*yFP)
    c = arctan(rho)
    decStar = arcsin(cos(c)*sin(decOpticalAxis)+(-xFP*sin(c)*cos(decOpticalAxis))/rho)
    raStar = raOpticalAxis + arctan2(yFP*sin(c), rho*cos(decOpticalAxis)*cos(c)+xFP*sin(decOpticalAxis)*sin(c))
    

    # Return the equatorial coordinates

    return raStar, decStar








def pixelToFocalPlaneCoordinates(xCCDpixel, yCCDpixel, pixelSize, ccdZeroPointX, ccdZeroPointY, CCDangle):

    """
    PUROSE: Given the (real-valued) pixel coordinates of the star on the CCD, compute the (x,y)
            coordinates in the FP' reference system (not the FP system!).

    INPUT: xCCDpixel     : x-coordinate of the star on the CCD  [pixel]
           yCCDpixel     : y-coordinate of the star on the CCD  [pixel]
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








def drawCCDsInSky(raOpticalAxis, decOpticalAxis, focalPlaneAngle, plateScale, pixelSize, nominal=True):

    """
    PURPOSE: Project and plot the 4 CCDs of 1 camera on the sky

    INPUT: raOpticalAxis:   right ascension of the optical axis [rad]
           decOpticalAxis:  declination of the optical axis [rad]
           focalPlaneAngle: angle between the Y_FP axis and the Y'_FP axis: gamme_FP  [rad]
           plateScale:      [arcsec/micron]
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

        # Get the corner coordinates in the FP' plane

        cornersXmm, cornersYmm = computeCCDcornersInFocalPlane(ccdCode, pixelSize)

        # Convert the FP' coordinates to equatorial coordinates

        ra, dec = inverseGnomonicProjectionFocalPlaneToSky(cornersXmm, cornersYmm, raOpticalAxis, decOpticalAxis, focalPlaneAngle, plateScale, pixelSize)

        # Repeat the coordinates of the 1st corner, to plot a nice closed loop
        # Convert from radians to degrees

        ra  = append(ra, ra[0]) * 180 / pi
        dec = append(dec, dec[0]) * 180 / pi
        
        plt.plot(ra, dec, c=color[ccdCode])

        # Overplot the row closest to the readout register with a thicker line

        plt.plot([ra[0], ra[1]], [dec[0], dec[1]], c=color[ccdCode], linewidth=3)


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


    plt.draw()

    # That's it

    return







def drawSubfieldInFocalPlane(ccdCode, xCCD, yCCD, subfieldSizeX, subfieldSizeY, pixelSize):

    """
    PURPOSE: Draw a subfield in the focal plane.
    """

    # Compute the position of the subfield in pixel coordinates, for the current CCD, 
    # disregarding the physical extend of the CCD

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








def drawPixelInFocalPlane(ccdCode, xCCD, yCCD, pixelSize):

    """
    PURPOSE: Plot a pixel from a particular CCD in the focal plane. The actual position in millimeter
             is shown as a red dot, while the pixel itself is drwan as a rectangle with edge pixelSize.

    INPUTS:  ccdCode:   for nominal camera: either 'A', 'B', 'C', or 'D'
                        for fast camer: either 'AF', 'BF', 'CF', 'DF'
             xCCDpix:   x-coordinate (column number, zero-based) of the pixel on the CCD  [pix]
             yCCDpix:   y-coordinate (row number, zero-based) of the pixel on the CCD  [pix]

    OUTPUTS: None
    """

    # Compute the position of the star in pixel coordinates, for the current CCD, 
    # disregarding the physical extend of the CCD

    zeroPointXmm = CCD[ccdCode]["zeroPointXmm"]
    zeroPointYmm = CCD[ccdCode]["zeroPointYmm"]
    ccdAngle     = CCD[ccdCode]["angle"]

    xFPprime, yFPprime = pixelToFocalPlaneCoordinates(xCCD, yCCD, pixelSize, zeroPointXmm, zeroPointYmm, ccdAngle)

    # Get the current axis

    currentAxis = plt.gca()
    currentAxis.add_patch( patches.Rectangle( (xFPprime, yFPprime), pixelSize / 1000.0, pixelSize / 1000.0, fill=False) )

    plt.plot(xFPprime, yFPprime, 'ro')
    plt.draw()

    return









def getCCDandPixelCoordinates(raStar, decStar, raOpticalAxis, decOpticalAxis, focalPlaneAngle,  \
                              focalLength, plateScale, pixelSize, nominal=True):

    """
    PURPOSE: Given the equatorial coordinates of a star, find out on which CCD it falls ('A', 'B', ...)
             and compute the pixel coordinates of the star on this CCD.

    INPUT: raStar:          right ascension of the star [rad]
           decStar:         declination of the star [rad]      
           raOpticalAxis:   right ascension of the optical axis [rad]
           decOpticalAxis:  declination of the optical axis [rad]
           focalPlaneAngle: angle between the Y_FP axis and the Y'_FP axis: gamme_FP  [rad]
           focalLength:     focal length of the telescope [m]
           plateScale:      [arcsec/micron]
           pixelSize:       [micrometer]
           nominal:         True for the nominal camera configuration, False for the fast cameras

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

    xFPrad, yFPrad = skyToAngularFocalPlaneCoordinates(raStar, decStar, raOpticalAxis, decOpticalAxis, focalPlaneAngle, plateScale, pixelSize)
    xFPmm, yFPmm = angularToPlanarFocalPlaneCoordinates(xFPrad, yFPrad, focalLength)

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









def getSkyCoordinates(ccdCode, xCCDpix, yCCDpix, plateScale, pixelSize, raOpticalAxis, decOpticalAxis, focalPlaneAngle):

    """
    PURPOSE: Return the sky coordinates of a pixel on the CCD in equatorial coordinates.

    INPUTS:  ccdCode:         for nominal camera: either 'A', 'B', 'C', 'D'
                              for fast camera: either 'AF', 'BF', 'CF', 'DF'
             xCCDpix:         x-coordinate (column number, zero-based) of the pixel on the CCD  [pix]
             yCCDpix:         y-coordinate (row number, zero-based) of the pixel on the CCD  [pix]
             plateScale:      [arcsec/micron]
             pixelSize:       [micrometer]
             raOpticalAxis:   right ascension of the optical axis [rad]
             decOpticalAxis:  declination of the optical axis [rad]
             focalPlaneAngle: angle between the Y_FP axis and the Y'_FP axis: gamme_FP  [rad]
 
    OUTPUTS: raStar:          right ascension of the star [rad]
             decStar:         declination of the star [rad]
    """

    # Compute the position of the star in pixel coordinates, for the current CCD, 
    # disregarding the physical extend of the CCD

    zeroPointXmm = CCD[ccdCode]["zeroPointXmm"]
    zeroPointYmm = CCD[ccdCode]["zeroPointYmm"]
    ccdAngle     = CCD[ccdCode]["angle"]

    xFPprime, yFPprime = pixelToFocalPlaneCoordinates(xCCDpix, yCCDpix, pixelSize, zeroPointXmm, zeroPointYmm, ccdAngle)

    # Convert the FP' coordinates to equatorial coordinates

    raStar, decStar = inverseGnomonicProjectionFocalPlaneToSky(xFPprime, yFPprime, raOpticalAxis, decOpticalAxis, focalPlaneAngle, plateScale, pixelSize)

    return raStar, decStar






def calculateSubfieldAroundCoordinates(raStar, decStar, subfieldSizeX, subfieldSizeY, focalLength, plateScale, pixelSize, \
                                       raOpticalAxis, decOpticalAxis, focalPlaneAngle, nominal=True):

    """
    PURPOSE: Calculates the location of the subfield such that the star with coordinates (raStar, decStar)
             is centered in the subfield.

    INPUTS:  raStar:          right ascension [rad]
             decStar:         declination [rad]
             subfieldSizeX:   full width (# of columns) of the subfield [pix]
             subfieldSizeY:   full height (#of rows) of the subfield [pix]
             focalLength:     focal length of the telescope [m]
             plateScale:      [arcsec/micron]
             pixelSize:       [micrometer]
             raOpticalAxis:   right ascension of the optical axis [rad]
             decOpticalAxis:  declination of the optical axis [rad]
             focalPlaneAngle: angle between the Y_FP axis and the Y'_FP axis: gamme_FP  [rad]
             nominal:         True for the nominal camera configuration, False for the fast cameras

    OUTPUTS: ccdCode: "A", "B", "C" or "D" if nominal=True, "AF", "BF", "CF" or "DF" otherwise
             xCCDpix: x-coordinate of the star in pixels (i.e. column number)
             yCCDpix: y-coordinate of the star in pixels (i.e. row number)

    REMARKS: - If the coordinates do not fall on any CCD, an error message is shown, followed by an exit(1)
             - If the star is too close to the edge for the given subfield size, and error message is shown,
               followed by an exit(1)
    """

    # Find out on which CCD the star falls, and the corresponding pixel coordinates

    ccdCode, xCCDpix, yCCDpix = getCCDandPixelCoordinates(raStar, decStar, raOpticalAxis, decOpticalAxis, \
        focalPlaneAngle, focalLength, plateScale, pixelSize, nominal)

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












def setSubfieldAroundCoordinates(sim, raStar, decStar, subfieldSizeX, subfieldSizeY, focalLength, plateScale, pixelSize, \
                                 raOpticalAxis, decOpticalAxis, focalPlaneAngle, nominal=True):
    
    """
    Calculates the location of the sub-field such that it is centred on the star 
    with the given sky coordinates.  Depending on the CCD (in nomincal mode:
    "A", "B", "C", or "D"; in fast mode: "AF", "BF", "CF", or "DF"), the 
    configuration file for the given simulation are adapted.  These include the
    pre-defined CCD position, the dimensions of the CCD (and also of the sub-field,
    although this is not affected by the calculations), the sub-field zeropoint
    and the exposure time. 

    INPUTS:  sim:             simulation for which the configuration file is adapted
             raStar:          right ascension of the star [radians]
             decStar:         declination [radians]
             subfieldSizeX:   width (i.e. number of columns) of the subiield [pixels]
             subfieldSizeY:   height (i.e. number of rows) of the sub-field [pixels]
             focalLength:     focal length of the telescope [m]
             plateScale:      Plate scale. [arcsec/micron]
             pixelSize:       [micrometer]
             raOpticalAxis:   right ascension of the optical axis [radians]
             decOpticalAxis:  declination of the optical axis [radians]
             focalPlaneAngle: orientation angle of the focal plane [radians]
             nominal:         True for the nominal camera configuration, False for the fast cameras

    OUTPUT: True if the CCD code (i.e. the pre-defined CCD position) could be
            determined, False otherwise 

    REMARKS: - If the coordinates do not fall on any CCD, an error message is shown, followed by an exit(1)
             - If the star is too close to the edge for the given subfield size, and error message is shown,
               followed by an exit(1)
    """
    
    # Compute the position of the subfield.
    # xPix and yPix are the CCD coordinates of the star, given a 4510x4510 CCD [colNumber, rowNumber].

    ccdCode, xPix, yPix = calculateSubfieldAroundCoordinates(raStar, decStar, subfieldSizeX, subfieldSizeY, focalLength, plateScale, pixelSize, raOpticalAxis, decOpticalAxis, focalPlaneAngle, nominal)
    
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

    sim["SubField/ZeroPointRow"] = str(int(xPix - subfieldSizeX/2))
    sim["SubField/ZeroPointColumn"] = str(int(yPix - subfieldSizeY/2))
    sim["SubField/NumRows"] = str(subfieldSizeX)
    sim["SubField/NumColumns"] = str(subfieldSizeY)

    # Set the exposure time, depending on fast vs nominal cams

    if nominal:
        sim["ObservingParameters/ExposureTime"] = 23
    else:
        sim["ObservingParameters/ExposureTime"] = 2.3
    
    # That's it

    return True


