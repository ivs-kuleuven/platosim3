from numpy import *
import numpy as np
from referenceFrames import *








def drawCCDsInSky(raPlatform, decPlatform, tiltAngle, azimuthAngle, focalPlaneAngle, focalLength, pixelSize, normal=True):

    """
    PURPOSE: Project and plot the 4 CCDs of 1 camera on the sky

    INPUT: raPlatform:      right ascension of the platform pointing axis             [rad]
           decPlatform:     declination of the platform pointing axis                 [rad]
           tiltAngle:       tilt angle of the telescope w.r.t. platform z-axis        [rad]                   
           azimuthAngle:    azimuth angle of the telescope on the platform            [rad]
           focalPlaneAngle: angle between the Y_TL axis and the Y_FP axis: gamma_FP   [rad]
           focalLength:     focal length of the camera                                [mm]
           pixelSize:       pixel size                                                [micron]
           normal:         True for the normal camera configuration, False for the fast cameras

    OUTPUT: None

    TODO: - Does not work yet for the fast cams
          - Does not take distoration into account yet
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

        ra, dec = focalPlaneToSkyCoordinates(cornersXmm, cornersYmm, raPlatform, decPlatform,   \
                                             tiltAngle, azimuthAngle, focalPlaneAngle, focalLength)

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

    INPUT:    sim:     instance of simulation class (see simulation.py)
              raStar:  right ascension of the star [rad]
              decStar: declination of the star     [rad]
    
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

    pixelSize        = float(sim["CCD/PixelSize"])
    raPlatform       = np.radians(float(sim["ObservingParameters/RApointing"]))
    decPlatform      = np.radians(float(sim["ObservingParameters/DecPointing"]))
    azimuthTelescope = np.deg2rad(float(sim["Telescope/AzimuthAngle"]))
    tiltTelescope    = np.deg2rad(float(sim["Telescope/TiltAngle"]))
    focalPlaneAngle  = np.radians(float(sim["Camera/FocalPlaneOrientation"]))
    focalLength      = float(sim["Camera/FocalLength"]) * 1000.0  # [m] -> [mm]
    ccdZeroPointX    = float(sim["CCD/OriginOffsetX"])
    ccdZeroPointY    = float(sim["CCD/OriginOffsetY"])
    ccdAngle         = np.radians(float(sim["CCD/Orientation"]))

    xFPmm, yFPmm = skyToFocalPlaneCoordinates(raStar, decStar, raPlatform, decPlatform, tiltTelescope, azimuthTelescope, \
                                              focalPlaneAngle, focalLength)

    if includeFieldDistortion:
        xFPmm, yFPmm = undistortedToDistortedFocalPlaneCoordinates(xFPmm, yFPmm)

    ccdCode, xCCD, yCCD = getCCDandPixelCoordinates(raStar, decStar, raPlatform, decPlatform, tiltTelescope, azimuthTelescope,  \
                                                    focalPlaneAngle, focalLength, pixelSize, includeFieldDistortion, normal)

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

    xFPmm, yFPmm = pixelToFocalPlaneCoordinates(xCCD, yCCD, pixelSize, zeroPointXmm, zeroPointYmm, ccdAngle)   # [mm]

    # Get the current axis

    currentAxis = plt.gca()
    currentAxis.add_patch( patches.Rectangle( (xFPmm, yFPmm), pixelSize / 1000.0, pixelSize / 1000.0, fill=False) )

    plt.plot(xFPmm, yFPmm, 'ro')
    plt.draw()

    return
