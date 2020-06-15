from simfile import SimFile
from astropy.io import fits
import referenceFrames as rf
import math

def fromHDF5toFITS(inputFilename, exposure, outputFilename):
    
    """
    PURPOSE: Convert the given exposure from the given PlatoSim HDF5 output file to a
             FITS file with a header compring the sub-field information.
    
    INPUT:
        - input_filename: Filename of the PlatoSim HDF5 file.
        - param exposure: Number of the exposure to include in the output FITS file.
        - output_filename: Filename of the FITS file comprising the given exposure 
                           and a header that contains the sub-field information.
    """

    sfile = SimFile(inputFilename)
    
    image =  sfile.getImage(exposure)
    
    ccdCode = int(sfile.getInputParameter("CCD", "Position"))
    subfieldZeropointRow = int(sfile.getInputParameter("SubField", "ZeroPointRow"))           # Sub-field zeropoint row [pixels]
    subfieldZeropointColumn = int(sfile.getInputParameter("SubField", "ZeroPointColumn"))     # Sub-field zeropoint column [pixels]
    
    ccdOriginOffsetX = sfile.getInputParameter("CCDPositions", "OriginOffsetX")[ccdCode - 1]      # CCD origin offset x-coordinate [mm]
    ccdOriginOffsetY = sfile.getInputParameter("CCDPositions", "OriginOffsetY")[ccdCode - 1]      # CCD origin offset y-coordinate [mm]
    ccdOrientationAngle = sfile.getInputParameter("CCDPositions", "Orientation")[ccdCode - 1]      # CCD orientation angle [degrees]
    
    pixelSize = sfile.getInputParameter("CCD", "PixelSize")    # Pixel size [µm]
    
    focalLength = sfile.getInputParameter("Camera/FocalLength", "ConstantValue") * 1000.0   # Focal length [mm]

    numRows = sfile.getInputParameter("SubField", "NumRows")           # Number of rows in the sub-field
    numColumns = sfile.getInputParameter("SubField", "NumColumns")     # Number of columns in the sub-field
    
    header = fits.Header()

    
    header["SIMPLE"] = "T"
    
    # Dimensions of the sub-field
    
    header["NAXIS"] = (2, "Dimensionality of the sub-field")
    header["NAXIS1"] = (numColumns, "Number of columns in the sub-field")
    header["NAXIS2"] = (numRows, "Number of rows in the sub-field")
    
    # CCD code
    
    header["CCD_ID"] = (ccdCode, "CCD code")
    
    # Focal length (this is needed for the conversion to field angles)
    
    header["FOCALLEN"] = (focalLength, "Focal length [mm]")
    
    # Linear coordinate transformation from sub-field to focal-plane coordinates
    
    header["ctype1"] = ("LINEAR", "Linear coordinate transformation")
    header["ctype2"] = ("LINEAR", "Linear coordinate transformation")

    # Focal-plane coordinates are expressed in mm

    header["CUNIT1"] = ("MM", "Target unit in the column direction (mm)")
    header["CUNIT2"] = ("MM", "Target unit in the row direction (mm)")
    
    # CCD origin + corresponding focal-plane coordinates
    
    # For the CCD origin, we know the coordinates in the source coordinate sytem
    # (i.e. the sub-field reference frame), from the sub-field zeropoint, and
    # can calculate the coordinates in the target coordinate system (i.e. in the 
    # focal-plane reference frame).  Hence, we use the CCD origin to be the 
    # reference point of the coordinate transformation (from the sub-field 
    # reference frame to the focal-plane reference frame):
    #   - the coordinates of the reference point in the source reference frame
    #     (i.e. in the sub-field reference frame) need to go in the CRPIXi
    #     keywords;
    #   - the coordinates of the reference point in the target reference frame
    #     (i.e. in the focal-plane reference frame) need to go in the CRVALi
    #     keywords.

    header["CRPIX1"] = (-subfieldZeropointColumn, "Sub-field column of the CCD origin [pixels]")
    header["CRPIX2"] = (-subfieldZeropointRow, "Sub-field row of the CCD origin [pixels]")

    crval1, crval2 = rf.pixelToFocalPlaneCoordinates(0, 0, pixelSize, ccdOriginOffsetX, ccdOriginOffsetY, math.radians(ccdOrientationAngle))

    header["CRVAL1"] = (crval1, "FP x-coordinate of the CCD origin [mm]")
    header["CRVAL2"] = (crval2, "FP y-coordinate of the CCD origin [mm]")
    
    # Pixel size
    
    cdelt = pixelSize / 1000.0
    header["CDELT1"] = (cdelt, "Pixel size in the x-direction [micron]")
    header["CDELT2"] = (cdelt, "Pixel size in the y-direction [micron]")
    
    # Orientation angle of the CCD

    header["crota2"] = (ccdOrientationAngle, "CCD orientation angle [degrees]")
    
    ccdOrientationAngleRadians = math.radians(ccdOrientationAngle)
    header["cd1_1"] = (cdelt * math.cos(ccdOrientationAngleRadians), "Pixel size x cos(CCD orientation angle)")
    header["cd1_2"] = (-cdelt * math.sin(ccdOrientationAngleRadians), "-Pixel size x sin(CCD orientation angle)")
    header["cd2_1"] = (cdelt * math.sin(ccdOrientationAngleRadians), "Pixel size x sin(CCD orientation angle)")
    header["cd2_2"] = (cdelt * math.cos(ccdOrientationAngleRadians), "Pixel size x cos(CCD orientation angle)")
    
    # Additional keywords for the future
    
    #header["INSTRUME"] or header["TELESCOP"] = ("", "Camera name")
    #header["SITENAME"] = ("", "Name of the test site")
    #header["EXPOSURE"] = ("", "Exposure time [s]")
    #header["DATE-LOC"] = ("", "Local time of observation")
    
    hdu = fits.PrimaryHDU(image)
    hdu.header = header
    
    hdu.writeto(outputFilename)
    