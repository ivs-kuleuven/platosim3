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

    print("Determining PSF sigma.")

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

    print("Sigma PSF = {0}".format(sigmaPSF))

    
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
            rowStar = int(rowPix[starNr])
            colStar = int(colPix[starNr])

            for row in range(rowStar-1, rowStar+2):
                
                if (row < 0) or (row > image.shape[0]-1): continue

                for col in range(colStar-1, colStar+2):

                    if (col < 0) or (col > image.shape[1]-1): continue
                    
                    weight = exp(-(pow(row-rowPix[starNr],2) + pow(col-colPix[starNr],2))/2.0/sigmaPSF/sigmaPSF)
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


