from numpy import *
from numpy.random import uniform
import os
import sys
import h5py



def computePSFsigma(psf, Nsubpixels, Nsamples=10000):

    """
    PURPOSE: Approximate the PSF with a 2D symmetric Gaussian distribution, to compute 
             its standard deviation, assuming that the barycenter is in the middle of the image. 
             This standard deviation will be used for the weighted mask photometry

    INPUT: psf:        2D numpy image containing the rotated PSF at subpixel level
           Nsubpixels: the number of subpixels per pixel (1D)
           Nsamples:   the number of monte carlo samples of the PSF to fit a Gaussian

    OUTPUT: sigma: the standard deviation of the symmetric PSF  [pix]

    """

    sumPSF = sum(psf)
    randomNumber = uniform(0.0, 1.0, Nsamples)

    rows = []
    cols = []

    for n in range(Nsamples):

        # Go through the PSF pixels <Nsamples> time, storing the pixel coordinates 
        # with higher PSF values more often than those with low values. 

        cumulative = 0.0;
        skipRestOfPSF = False

        # Sum up the PSF values until you exceed the generated random number
        # Then keep the (row,col) value that made the sum exceed the threshold,
        # and skip the rest of the PSF pixels.

        for i in range(psf.shape[0]):
            if skipRestOfPSF is True: break
            for j in range(psf.shape[1]):
                if skipRestOfPSF is True: break
                cumulative = cumulative + psf[i,j];
                if cumulative >= randomNumber[n] * sumPSF:
                    rows.append(i)
                    cols.append(i)
                    skipRestOfPSF = True
                    break

    # Compute the standard deviation of the PSF as the average of the sigma in the row direction
    # and the sigma in the colum  direction    
    
    sigmaPSF = (array(rows).std() + array(cols).std())/2.0                       # [supixels]
    
    # Convert the stdev of the PSF from subpixels to pixels

    sigmaPSF = sigmaPSF / Nsubpixels                                             # [pixels]

    # That's it!

    return sigmaPSF








def computePSFBarycenterOffset(psf, Npixels):

    """
    PURPOSE: Compute the barycenter coordinates of the given PSF. The zero-point is the center of the psf.

    INPUT: psf:     2D numpy image (Nrows x Ncolumns) containing the rotated PSF at subpixel level
           Npixels: the number of pixels in 1 dimension of the square PSF (e.g. 8)

    OUTPUT: (deltaRow, deltaColumn): the pixel coordinates of the barycenter of the PSF
                                     relative to the center of the pixel image   [pix]
    """

    Nrows, Ncols = psf.shape

    rowIndices = arange(Nrows)
    colIndices = arange(Ncols) 

    baryRow = sum(psf * rowIndices) / psf.sum()
    baryCol = sum(psf * colIndices.reshape(Ncols,1)) / psf.sum()

    Nsubpixels = Nrows / Npixels
    deltaRow = (baryRow - Nrows / 2.0) / Nsubpixels
    deltaCol = (baryCol - Ncols / 2.0) / Nsubpixels

    return (deltaRow, deltaCol)










def photometry(inputFilePath, outputFilePath, sigmaPSF, verbose=False):

    """
    PURPOSE: Given a platosim simulation output file, this function loops over all
             stored images, and performs a weighted aperture photometry 

    INPUT: inputFilePath:  path of the input hdf5 file
           outputFilePath: path of the output hdf5 file
           sigmaPSF: the standard deviation of the PSF (assumed to be symmetrical) [pix]
           verbose: if False, do not print info messages on screen

    OUTPUT: None

    """


    # Open the input and the output HDF5 files

    inputFile = h5py.File(inputFilePath, "r")
    
    if not os.path.isfile(outputFilePath):
        outputFile = h5py.File(outputFilePath, "w")
    else:
        raise IOError("Output file already exists")


    # Create the necessary groups in the output hdf5 file

    photometryGroup = outputFile.create_group("/Photometry");

    # Copy the time points of the input HDF5 file to the output file.

    time = array(inputFile["/StarPositions/Time"])
    photometryGroup.create_dataset("time", data=time)

    # Compute the offset of the barycenter of the PSF

    psf = array(inputFile["/PSF/rotatedPSF"])                                           # [subpix]
    if inputFile["/InputParameters/PSF"].attrs["Model"] == "FromFile":
        Npixels = int(inputFile["/InputParameters/PSF/FromFile"].attrs["NumberOfPixels"])
    else:
        Npixels = int(inputFile["/InputParameters/PSF/Gaussian"].attrs["NumberOfPixels"])

    deltaRow, deltaCol = computePSFBarycenterOffset(psf, Npixels)                       # [pix]


    # Collect the relevant input parameters from the HDF5 file

    gain = inputFile["/InputParameters/CCD/"].attrs["Gain"]                             # [e-/ADU]
    quantumEfficiency = inputFile["/InputParameters/CCD"].attrs["QuantumEfficiency"]    # [e-/phot]
    sigmaRON = inputFile["/InputParameters/CCD"].attrs["ReadoutNoise"]                  # [e-/pix]
    

    # Extract the sky background. Convert from [phot/pix/exposure] to [e-/pix/exposure]
    # Electrons are always integer numbers, so round down.

    skyBackground = array(inputFile["/Background/skyBackground"])                       # [phot/pix/exposure]
    skyBackground = floor(skyBackground * quantumEfficiency)                            # [e-/pix/exposure]


    # Get the magnitudes of all the stars found in any of the images. The number of stars 
    # can differ from image to image because some stars may jitter out of the subfield. 
    # Extract arrays from the HDF5 file, but put it in a more convenient dictionary 'Vmag'.
    # E.g. magnitude[12445] contains the magnitude of the star with ID 12445.

    id = array(inputFile["StarCatalog/starIDs"])
    mag = array(inputFile["StarCatalog/Vmag"])
    magnitude = dict([(id[n], mag[n]) for n in range(len(mag))])

    # Loop over all exposure, and apply weighted mask photometry on each image

    Nexposures = inputFile["/InputParameters/ObservingParameters/"].attrs["NumExposures"];

    if verbose:
        print("Looping over all images in HDF5 file.")

    for imageNr in range(Nexposures):

        if verbose:
            print("Image # {0}".format(imageNr))

        # Read the bias and smearing map, and the image itself

        image = array(inputFile["/Images/image{0:06d}".format(imageNr)])
        smearingMap = array(inputFile["/SmearingMaps/smearingMap{0:06d}".format(imageNr)])
        biasMap = array(inputFile["/BiasMaps/biasMap{0:06d}".format(imageNr)])


        # Estimate the bias [ADU] and subtract it from the image and the smearing map

        bias = floor(biasMap.mean())
        image -= bias
        smearingMap -= bias

        if verbose:
            print("    Subtracted bias level of {0} ADU".format(bias))

        # Correct for open shutter smearing using the smearing maps
        # meanSmearing contains a smearing value for each column

        meanSmearing = smearingMap.mean(axis=0)
        image -= meanSmearing

        if verbose:
            print("    Corrected for open-shutter smearing")

        # Convert from [ADU] to [electrons] using the gain

        image *= gain

        if verbose:
            print("    Converted from [ADU] to [electrons] using a Gain of {0} e-/ADU".format(gain))

        # Correct for geometrical vignetting
        # TODO


        # Correct for the flatfield

        flatfield = array(inputFile["Flatfield/PRNU"])
        image /= flatfield;

        if verbose:
            print("    Corrected for PRNU")
        
        # Loop over all stars in this image, and do weighted aperture photometry

        exposureGroupName = "/Photometry/Exposure{0:06d}".format(imageNr)
        exposureGroup = outputFile.create_group(exposureGroupName)

        starID    = array(inputFile["StarPositions/Exposure{0:06d}/starID".format(imageNr)])
        colPix    = array(inputFile["StarPositions/Exposure{0:06d}/colPix".format(imageNr)])
        rowPix    = array(inputFile["StarPositions/Exposure{0:06d}/rowPix".format(imageNr)])
        inputFlux = array(inputFile["StarPositions/Exposure{0:06d}/flux".format(imageNr)]) 

        # The input flux stored in the HDF5 file, are expressed in [phot/exposure],
        # and the passband and the transmission efficiency of the telescope are already included
        # in the value. We still need to convert from [photons/exposure] to [e-/exposure].

        inputFlux = inputFlux * quantumEfficiency
        

        # Loop over all stars in the image, and extract the flux

        maskSize = zeros(len(starID))
        estimatedFlux = zeros(len(starID))
        varEstimatedFlux = zeros(len(starID))
        Vmag = zeros(len(starID))

        if verbose:
            print ("    Looping over all stars in image to do weighted mask photometry")
            print ("        Using background level of {0} e-/pix/exposure".format(skyBackground[imageNr]))

        for starNr in range(len(starID)):
            
            Vmag[starNr] = magnitude[starID[starNr]]
            estimatedFlux[starNr] = 0.0
            rowStar = int(rowPix[starNr] + deltaRow)         # deltaRow, deltaCol are the PSF barycenter offsets
            colStar = int(colPix[starNr] + deltaCol)

            for row in range(rowStar-1, rowStar+2):
                
                if (row < 0) or (row > image.shape[0]-1): continue

                for col in range(colStar-1, colStar+2):

                    if (col < 0) or (col > image.shape[1]-1): continue
                    
                    weight = exp(-(pow(row-rowPix[starNr]-deltaRow,2) + pow(col-colPix[starNr]-deltaCol,2))/2.0/sigmaPSF/sigmaPSF)
                    estimatedFlux[starNr] += weight * (image[row, col] - skyBackground[imageNr])
                    varEstimatedFlux[starNr] += weight*weight * (image[row, col] + sigmaRON * sigmaRON)
                    maskSize[starNr] +=1 


        # Compute the S/N for each star

        snr = estimatedFlux / sqrt(varEstimatedFlux)

        # Save all computed arrays for this image to the HDF5 file

        exposureGroup.create_dataset("starID",           data=starID)
        exposureGroup.create_dataset("estimatedFlux",    data=estimatedFlux)     # [e-/exposure]
        exposureGroup.create_dataset("varEstimatedFlux", data=varEstimatedFlux)  # [e-/exposure]
        exposureGroup.create_dataset("inputFlux",        data=inputFlux)         # [e-/exposure]
        exposureGroup.create_dataset("SNR",              data=snr)            
        exposureGroup.create_dataset("Vmag",             data=Vmag)                        
        exposureGroup.create_dataset("maskSize",         data=maskSize)


    # All images have now been processed. Close the HDF5 files.

    inputFile.close();
    outputFile.close();
    
    # That's it!

    #return estimatedFlux, varEstimatedFlux, snr, inputFlux, Vmag, maskSize






def getPhotometryTimeSeries(photometryFile, starID):

    """
     PURPOSE: extract the flux time series of star with a given identifier.

     INPUT: photometryFile: an HDF5 output file written by the photometry() function above
            starID:  star identifier (integer, e.g. 9789)

     OUTPUT: time: a numpy array containing the time points [s]
             flux: a numpy array containing the flux points [electrons/exposure]

     REMARK: To find out which star identifiers are in the photometry file, look in the HDF5 simulation
             output file of PlatoSim: 
             allStarIDs = array(platosimOutputFile["StarCatalog/starIDs"])
    """

    photFile = h5py.File(photometryFile)
    allTimePoints = array(photFile["/Photometry/time"])
    Nimages = len(allTimePoints)

    time = []
    flux = []

    for k in range(Nimages):
        allStarIDsInImage = array(photFile["/Photometry/Exposure{0:06d}/starID".format(k)])
        if starID in allStarIDsInImage:
            estimatedFlux = array(photFile["/Photometry/Exposure{0:06d}/estimatedFlux".format(k)])
            flux.append(estimatedFlux[where(allStarIDsInImage==starID)][0])
            time.append(array(photFile["/Photometry/time"])[k])

    flux = array(flux)
    time = array(time)

    photFile.close()

    return time, flux
    

