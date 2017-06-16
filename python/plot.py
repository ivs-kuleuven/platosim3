from numpy import *
import numpy as np
from referenceFrames import *











def drawCCDsInSkyMollweide(fig, raPlatform, decPlatform, tiltAngle, azimuthAngle, focalPlaneAngle, focalLength, pixelSize, normal=True):

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

        ra, dec = focalPlaneToSkyCoordinates(cornersXmm, cornersYmm, raPlatform, decPlatform,   \
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







def drawStarsInSkyMollweide(fig, ra, dec ):

    """
    PURPOSE: Project and plot the stars with the given right ascension and declination on the sky

    INPUT: ra:      right ascension of the stars             [rad]
           dec:     declination of the areA                 [rad]

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
        ccdCodes = ['1', '2', '3', '4']
    else:
        ccdCodes = ['1F', '2F', '3F', '4F']


    # Set up the colors to be used to draw each CCD. 
    # Different CCDs have different colors.

    color = {'1': 'b', '1F': 'b', '2': 'r', '2F': 'r', '3': 'g', '3F': 'g', '4': 'k', '4F': 'k'}


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

    distortionCoefficients = None
    if (sim["Camera/IncludeFieldDistortion"] == "yes")  or (sim["Camera/IncludeFieldDistortion"] == "1"):
        includeFieldDistortion = True
        FIELD_DISTORTION["Coeff"] = sim["Camera/FieldDistortion/Coefficients"]
        FIELD_DISTORTION["InverseCoeff"] = sim["Camera/FieldDistortion/InverseCoefficients"]
        distortionCoefficients = sim["Camera/FieldDistortion/Coefficients"]
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

    #ccdCode, xCCD, yCCD = getCCDandPixelCoordinates(raStar, decStar, raPlatform, decPlatform, tiltTelescope, azimuthTelescope,  \
    #                                                focalPlaneAngle, focalLength, pixelSize, includeFieldDistortion, FIELD_DISTORTION["Coeff"], normal)
    ccdCode, xCCD, yCCD = getCCDandPixelCoordinates(raStar, decStar, raPlatform, decPlatform, tiltTelescope, azimuthTelescope, focalPlaneAngle, focalLength, pixelSize, includeFieldDistortion, distortionCoefficients, normal)

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
