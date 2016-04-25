from numpy import *
from numpy.random import uniform
import os
import sys
import h5py









def photometry(inputFilePath, outputFilePath):

    """
    PURPOSE: Given a platosim simulation output file, this function loops over all
             stored images, and performs a weighted aperture photometry 

    INPUT: inputFilePath:  path of the input hdf5 file
           outputFilePath: path of the output hdf5 file

    OUTPUT:

    """


    # Open the input and the output HDF5 files

    inputFile = h5py.File(inputFilePath, "r")
    
    if not os.path.isfile(outputFilePath):
        outputFile = h5py.File(outputFilePath, "w")
    else:
        raise IOError("Output file already exists")


    # Create the necessary groups in the output hdf5 file

    outputFile.create_group("/Photometry");

    # Read the PSF image (at subpixel level)

    psf = array(inputFile["/PSF/rebinnedPSFsubPixel"])

    # Approximate the PSF with a 2D symmetric gaussian distribution, to compute its standard deviation,
    # assuming that the barycenter is in the middle of the image. This standard deviation will be used
    # for the weighted mask photometry

    Nsamples = 10000
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

    Nsubpixels = inputFile["/InputParameters/SubField/"].attrs["SubPixels"]
    sigmaPSF = sigmaPSF / Nsubpixels                                             # [pixels]

    print("Sigma PSF: {0}".format(sigmaPSF))

    
    # Collect info from the HDF5 file necessary for the computations, but independent of the image Nr.

    gain = inputFile["/InputParameters/CCD/"].attrs["Gain"]
    quantumEfficiency = inputFile["/InputParameters/CCD"].attrs["QuantumEfficiency"]
    skyBackground = array(inputFile["/Background/skyBackground"])                       # [phot/pix/exposure]


    # Loop over all exposure, and apply weighted mask photometry on each image

    Nexposures = inputFile["/InputParameters/ObservingParameters/"].attrs["NumExposures"];

    for imageNr in range(Nexposures):

        # Read the bias and smearing map, and the image itself

        image = array(inputFile["/Images/image{0:06d}".format(imageNr)])
        smearingMap = array(inputFile["/SmearingMaps/smearingMap{0:06d}".format(imageNr)])
        biasMap = array(inputFile["/BiasMaps/biasMap{0:06d}".format(imageNr)])


        # Estimate the bias [ADU] and subtract it from the image and the smearing map

        bias = biasMap.mean()
        image -= bias
        smearingMap -= bias

        print("Bias level of image #{0}: {1}".format(imageNr, bias))

        # Correct for open shutter smearing using the smearing maps
        # meanSmearing contains a smearing value for each column

        meanSmearing = smearingMap.mean(axis=0)
        image -= meanSmearing

        # Convert from [ADU] to [electrons] using the gain

        image *= gain

        print("Gain: {0}".format(gain))

        # Convert from [electrons] to [photons] using the QE

        image /= quantumEfficiency;

        print("QE: {0}".format(quantumEfficiency))

        # Correct for geometrical vignetting
        # TODO


        # Correct for the flatfield

        flatfield = array(inputFile["Flatfield/PRNU"])
        image /= flatfield;

        
        # Correct for the background (stored in the HDF5 files in [photons/pix/exposure])
        # The background value stored in the input HDF5 file is the flux before entering the telescope.

        image -= skyBackground[imageNr];


        # Loop over all stars in this image, and do weighted aperture photometry

        exposureGroupName = "/Photometry/Exposure{0:06d}".format(imageNr)
        exposureGroup = outputFile.create_group(exposureGroupName)

        colPix    = array(inputFile["StarPositions/Exposure{0:06d}/colPix".format(imageNr)])
        rowPix    = array(inputFile["StarPositions/Exposure{0:06d}/rowPix".format(imageNr)])
        inputFlux = array(inputFile["StarPositions/Exposure{0:06d}/flux".format(imageNr)])
        starID    = array(inputFile["StarPositions/Exposure{0:06d}/starID".format(imageNr)])

        maskSize = zeros(len(starID))
        estimatedFlux = zeros(len(starID))

        for starNr in range(len(starID)):
            
            estimatedFlux[starNr] = 0.0
            rowStar = int(rowPix[starNr])
            colStar = int(colPix[starNr]

            for row in range(rowStar-1, rowStar+2):
                
                if (row < 0) or (row > image.shape[0]-1): continue

                for col in range(colStar-1, colStar+2):

                    if (col < 0) or (col > image.shape[1]-1): continue
                    
                    weight = exp(-(pow(row-rowPix[starNr],2) + pow(col-colPix[starNr],2))/2.0/sigmaPSF/sigmaPSF)
                    estimatedFlux[starNr] += weight * image[row, col]
                    maskSize[starNr] +=1 

            estimatedFlux[starNr] *= normalizationFactorNumerator / normalizationFactorDenominator


        exposureGroup.create_dataset("estimatedFlux", data=estimatedFlux)
        exposureGroup.create_dataset("inputFlux", data=inputFlux)
        exposureGroup.create_dataset("maskSize", data=maskSize)


    inputFile.close();
    outputFile.close();
    
    return estimatedFlux, inputFlux, maskSize


