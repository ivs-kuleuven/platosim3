#include "ClosedLoopDetectorClasses.h"


/**
 * \brief override the takeExposure function to include the receiving of windowPosition from server 
 *        
 */
double ClosedLoopDetectorWithMappedPSF::takeExposure(int exposureNr, double startTime, double exposureTime)
{
    // Advance the internal clock until the given start time

    internalTime = startTime;

    // check if there are new updates for the window position

    if (getWindowPositionFromServer)
    {
        // receive window postion
        std::tuple<bool, uint, uint, uint, uint, double> windowPositionTuple = getNewWindowPosition(exposureTime);

        // set the window position according to the new message
        setNewWindowPosition(windowPositionTuple);
    }

    // Integration of point sources and background, taking into account jitter + drift.

    Log.info("ClosedLoopDetectorWithMappedPSF: Integrating light for exposure " + to_string(exposureNr) + " with exposure time = " + to_string(exposureTime));

    integrateLight(exposureNr, startTime, exposureTime);

    // Include noise effects like readout noise, photon noise, full well saturation, etc.
    // Note: readOut() needs the exposure time to compute the open shutter smearing.

    Log.info("ClosedLoopDetectorWithMappedPSF: Adding noise effects to exposure " + to_string(exposureNr));

    readOut(exposureTime);

    // Write the CCD subfield, the bias map, and the smearing map to the HDF5 file

    Log.debug("ClosedLoopDetectorWithMappedPSF: Writing PixelMap, smearing map, and bias map #" + to_string(exposureNr) + " to HDF5 file.");

    writePixelMapsToHDF5(exposureNr);

    // Advance the internal clock

    internalTime += exposureTime + readoutTimeBeforeNextExposure;

    return internalTime;
}


/**
 * \brief override the writePixelMapsToHDF5 function to include the sending of the imagettes to client 
 *        
 */
void ClosedLoopDetectorWithMappedPSF::writePixelMapsToHDF5(int exposureNr)
{
    // Compose the image name

    stringstream myStream;
    myStream << "image" << setfill('0') << setw(6) << exposureNr;
    string imageName = myStream.str();

    // Add the image to the "Images" group

    if (!includeQuantisation)
    {
        // Write the float array to HDF5

        hdf5File.writeArray("/Images", imageName, pixelMap);
    }
    else
    {
        // Write the pixel maps as 2-byte (16 bit) unsigned short integers.
        // As a safety check, first check that the extrema of the map are indeed
        // within the boundaries of such a data type.
       
        if((pixelMap.min() < 0) || (pixelMap.max() >= (1 << 16)))
        {
            throw ConfigurationException("ClosedLoopDetectorWithMappedPSF: quantisation was applied but pixel map values are not in [0, 2^16[");
        }

        // Convert the float matrix to an unsigned uint16_t matrix

        arma::Mat<uint16_t> uintMap = arma::conv_to<arma::Mat<uint16_t>>::from(pixelMap);
        hdf5File.writeArray("/Images", imageName, uintMap);
    }

    // send the imagette to the client connected via tcp connection

    if (sendImagettesToClient)
    {
        Log.info("ClosedLoopDetectorWithMappedPSF: attempt to send the imagette to the client");

        sendImagetteToClient(&pixelMap, exposureNr);

    }


    if (numRowsSmearingMap != 0)
    {
        // Clear the string stream and compose the smearing map name

        myStream.str(string());      // insert empty string
        myStream.clear();            // clear eof bit

        myStream << "smearingMap" << setfill('0') << setw(6) << exposureNr;
        string smearingMapName = myStream.str();


        // Add the smearing map to the "SmearingMaps" group

        if (!includeQuantisation)
        {
            // Write the float array to HDF5

            //hdf5File.writeArray("/SmearingMaps", smearingMapName, smearingMap);
        }
        else
        {
            if ((smearingMap.min() < 0) || (smearingMap.max() >= (1 << 16)))
            {
                throw ConfigurationException("Detector: quantisation was applied but smearing map values are not in [0, 2^16[");
            }

            // Convert the float matrix to an unsigned uint16_t matrix

            arma::Mat<uint16_t> uintMap = arma::conv_to<arma::Mat<uint16_t>>::from(smearingMap);
            //hdf5File.writeArray("/SmearingMaps", smearingMapName, uintMap);
        }
        
    }

   // Clear the string stream and compose the bias map name

    myStream.str(string());      // insert empty string
    myStream.clear();            // clear eof bit

    myStream << "biasMap" << setfill('0') << setw(6) << exposureNr;
    string biasMapName = myStream.str();

    // Add the bias map to the "BiasMaps" group

    if (!includeQuantisation)
    {
        // Write the float array to HDF5

        //hdf5File.writeArray("/BiasMapsLeft", biasMapName, biasMapLeft);
        //hdf5File.writeArray("/BiasMapsRight", biasMapName, biasMapRight);
    }
    else
    {
        if ((biasMapLeft.min() < 0) || (biasMapLeft.max() >= (1 << 16)))
        {
            throw ConfigurationException("ClosedLoopDetectorWithMappedPSF: quantisation was applied but pixel values in the left bias map are not in [0, 2^16[");
        }

        if ((biasMapRight.min() < 0) || (biasMapRight.max() >= (1 << 16)))
        {
            throw ConfigurationException("ClosedLoopDetectorWithMappedPSF: quantisation was applied but pixel values in the right bias map are not in [0,2^16[");
        }

        // Convert the float matrix to an unsigned uint16_t matrix

        arma::Mat<uint16_t> uintMap = arma::conv_to<arma::Mat<uint16_t>>::from(biasMapLeft);
        //hdf5File.writeArray("/BiasMapsLeft", biasMapName, uintMap);

        uintMap = arma::conv_to<arma::Mat<uint16_t>>::from(biasMapRight);
        //hdf5File.writeArray("/BiasMapsRight", biasMapName, uintMap);
    }
    

    // Clear the string stream and compose the throughput map name

    myStream.str(string());      // insert empty string
    myStream.clear();            // clear eof bit

    myStream << "throughputMap" << setfill('0') << setw(6) << exposureNr;
    string throughputMapName = myStream.str();

    // Add the throughput map to the "ThroughputMaps" group

    //hdf5File.writeArray("/ThroughputMaps", throughputMapName, throughputMap);
}

 /**
 * \brief: variant of the generateFlatfieldMap without the writing of the maps to the hdf5 file - this is temporarity
 */

void ClosedLoopDetectorWithMappedPSF::generateFlatfieldMap()
{

    Log.info("Detector: generating flatfield map.");

    // Random number generation

    mt19937 flatfieldGenerator(flatfieldSeed);
    normal_distribution<double> flatfieldDistribution(0.0, 1.0);

    // Double the dimensions (this is necessary because of the behaviour of the Fourier transforms)
    // (this is a bit inconvenient as we are working at sub-pixel level -> to be investigated)

    int Nrows = 2 * numRowsPixelMap * numSubPixelsPerPixel;
    int Ncolumns = 2 * numColumnsPixelMap * numSubPixelsPerPixel;

    arma::cx_fmat evenMap = arma::cx_fmat(Nrows, Ncolumns);

    for(unsigned int row = 0; row < Nrows; row++)
    {
        for(unsigned int column = 0; column < Ncolumns; column++)
        {
            // Fourier space: generate white noise and include 1/f dependency
            // (Note: see https://en.wikipedia.org/wiki/Pink_noise#Generalization_to_more_than_one_dimension)

            evenMap(row, column) = flatfieldDistribution(flatfieldGenerator) / (pow(row, 2) + std::pow(column, 2) + 1);
        }
    }

    // Take the real part of the inverse Fourier transform

    evenMap = arma::ifft2(evenMap);
    arma::fmat realMap = arma::real(evenMap);

    // Cut out the appropriate part

    unsigned int numRowsFlatfield = Nrows / 2;
    unsigned int numColumnsFlatfield = Ncolumns / 2;
    
    flatfieldMap(arma::span::all, arma::span::all) = realMap(arma::span(0, numRowsFlatfield - 1), arma::span(0, numColumnsFlatfield - 1));
    flatfieldMap.reshape(numRowsFlatfield * numColumnsFlatfield, 1);

    // Normalisation
    //  - divide by mean and subtract 1.0 -> mean = 0.0
    //  - scale such that std.dev. = flatfield RMS and mean = 0.0
    //  - add 1.0

    flatfieldMap /= arma::mean(flatfieldMap.col(0));
    flatfieldMap -= 1;
    double scale = flatfieldNoiseRMS / arma::stddev(flatfieldMap.col(0));
    flatfieldMap *= scale;
    flatfieldMap += 1;

    flatfieldMap.reshape(numRowsFlatfield, numColumnsFlatfield);

    // Write the result to the HDF5 output file

    //hdf5File.writeArray("/Flatfield", "IRNU", flatfieldMap);

    // Rebin the intra-pixel flatfield to the pixel flatfield (IRNU -> PRNU)
    // and also write this array to the HDF5 outputfile. This PRNU array is not used
    // in the remainder of the simulation.

    arma::Mat<float> prnu(numRowsPixelMap, numColumnsPixelMap, arma::fill::zeros);

    for (unsigned int row = 0; row < numRowsPixelMap; row++)
    {
        for (unsigned int column = 0; column < numColumnsPixelMap; column++)
        {
            const unsigned int beginRow = row * numSubPixelsPerPixel;
            const unsigned int beginCol = column * numSubPixelsPerPixel;
            const unsigned int endRow = (row + 1) * numSubPixelsPerPixel - 1;
            const unsigned int endCol = (column + 1) * numSubPixelsPerPixel - 1;

            prnu(row, column) = arma::accu(flatfieldMap.submat(beginRow, beginCol, endRow, endCol))
                                / (numSubPixelsPerPixel * numSubPixelsPerPixel);
        }
    }

    // Write the result to the HDF5 output file

    Log.debug("Detector: writing PRNU to HDF5");

    //hdf5File.writeArray("/Flatfield", "PRNU", prnu);
}


/**
 * \brief change the window position and adapt all relevant pixel maps, if there was a new message from the server
 *        
 */
void ClosedLoopDetectorWithMappedPSF::setNewWindowPosition(std::tuple<bool, uint, uint, uint, uint, double> windowPositionTuple)
{
    // check whether a position change is necessary
    if (get<0>(windowPositionTuple))
    {
        // change the position of the subfield

        numRowsPixelMap         = get<1>(windowPositionTuple);
        numColumnsPixelMap      = get<2>(windowPositionTuple);

        subFieldZeroPointRow    = get<3>(windowPositionTuple);
        subFieldZeroPointColumn = get<4>(windowPositionTuple);
        
        orientationAngle      = deg2rad(get<5>(windowPositionTuple));
        
        Log.info("ClosedLoopDetectorWithMappedPSF: Changed numRowsPixelMap to: " + to_string(numRowsPixelMap));
        Log.info("ClosedLoopDetectorWithMappedPSF: Changed numColumnsPixelMap to: " + to_string(numColumnsPixelMap));

        Log.info("ClosedLoopDetectorWithMappedPSF: Changed subFieldZeroPointColumn to: " + to_string(subFieldZeroPointColumn));
        Log.info("ClosedLoopDetectorWithMappedPSF: Changed subFieldZeroPointRow to: " + to_string(subFieldZeroPointRow));

        Log.info("ClosedLoopDetectorWithMappedPSF: Changed orientationAngle to: " + to_string(orientationAngle));
            

        pixelMap.set_size(numRowsPixelMap, numColumnsPixelMap);

        biasMapLeft.set_size(numRowsBiasMap, numColumnsBiasMap);

        biasMapRight.set_size(numRowsBiasMap, numColumnsBiasMap);

        smearingMap.set_size(numRowsSmearingMap, numColumnsPixelMap);

        throughputMap.set_size(numRowsPixelMap, numColumnsPixelMap);

        throughputMap.ones(numRowsPixelMap, numColumnsPixelMap);

        // If we are going to apply open-shutter smearing, we have to know which pixels are within
        // the FOV (relevant only in case of mechanical vignetting).  When mechanical vignetting is
        // disabled, all pixels of the detector are inside the FOV.
            
        if(includeOpenShutterSmearing)
        {
            // Mechanical vignetting map:
            //  - no mechanical vignetting: all pixels of the sub-field inside FOV -> all values set to one
            //  - mechanical vignetting: set value of the pixels in the sub-field outside FOV to zero (others should be one) -> on creation of the throughput map
            
            mechanicalVignettingMask.ones(numRowsPixelMap, numColumnsPixelMap);

            // Number of exposed rows in each column:
            // - no mechanical vignetting: all exposed rows inside FOV (numRows - firstRowExposed)
            // - mechanical vignetting: count the exposed rows (i.e. from firstRowExposed) that are inside FOV
            
            numExposedRowsInFOV.zeros(numColumnsPixelMap);

            if(!includeMechanicalVignetting)
            {
                numExposedRowsInFOV.fill(numRows - firstRowExposed);
            }

            numRowsSubPixelMap = numRowsPixelMap * numSubPixelsPerPixel; // TODO Add edge pixels
            numColumnsSubPixelMap = numColumnsPixelMap * numSubPixelsPerPixel; // TODO Add edge pixels
            
            // Allocate memory for the different maps

            subPixelMap.zeros(numRowsSubPixelMap, numColumnsSubPixelMap);

            flatfieldMap.ones(numRowsSubPixelMap, numColumnsSubPixelMap);

            if(includeFlatfield)
            {
                // Generate the flatfield map

                generateFlatfieldMap();
            }

            // set the psf again for the new
            setPsfForSubfield();

        }
    }
}



// ------------------------------------------------------------------------------------------------------------- //




/**
 * \brief override the takeExposure function to include the receiving of windowPosition from server 
 *        
 */
double ClosedLoopDetectorWithAnalyticNonGaussianPSF::takeExposure(int exposureNr, double startTime, double exposureTime)
{
    // Advance the internal clock until the given start time

    internalTime = startTime;

    // check if there are new updates for the window position

    if (getWindowPositionFromServer)
    {
        // receive window postion
        std::tuple<bool, uint, uint, uint, uint, double> windowPositionTuple = getNewWindowPosition(exposureTime);

        // set the window position according to the new message
        setNewWindowPosition(windowPositionTuple);
    }

    // Integration of point sources and background, taking into account jitter + drift.

    Log.info("ClosedLoopDetectorWithMappedPSF: Integrating light for exposure " + to_string(exposureNr) + " with exposure time = " + to_string(exposureTime));

    integrateLight(exposureNr, startTime, exposureTime);

    // Include noise effects like readout noise, photon noise, full well saturation, etc.
    // Note: readOut() needs the exposure time to compute the open shutter smearing.

    Log.info("ClosedLoopDetectorWithMappedPSF: Adding noise effects to exposure " + to_string(exposureNr));

    readOut(exposureTime);

    // Write the CCD subfield, the bias map, and the smearing map to the HDF5 file

    Log.debug("ClosedLoopDetectorWithMappedPSF: Writing PixelMap, smearing map, and bias map #" + to_string(exposureNr) + " to HDF5 file.");

    writePixelMapsToHDF5(exposureNr);

    // Advance the internal clock

    internalTime += exposureTime + readoutTimeBeforeNextExposure;

    return internalTime;
}


/**
 * \brief override the writePixelMapsToHDF5 function to include the sending of the imagettes to client 
 *        
 */
void ClosedLoopDetectorWithAnalyticNonGaussianPSF::writePixelMapsToHDF5(int exposureNr)
{
    // Compose the image name

    stringstream myStream;
    myStream << "image" << setfill('0') << setw(6) << exposureNr;
    string imageName = myStream.str();

    // Add the image to the "Images" group

    if (!includeQuantisation)
    {
        // Write the float array to HDF5

        hdf5File.writeArray("/Images", imageName, pixelMap);
    }
    else
    {
        // Write the pixel maps as 2-byte (16 bit) unsigned short integers.
        // As a safety check, first check that the extrema of the map are indeed
        // within the boundaries of such a data type.
       
        if((pixelMap.min() < 0) || (pixelMap.max() >= (1 << 16)))
        {
            throw ConfigurationException("ClosedLoopDetectorWithMappedPSF: quantisation was applied but pixel map values are not in [0, 2^16[");
        }

        // Convert the float matrix to an unsigned uint16_t matrix

        arma::Mat<uint16_t> uintMap = arma::conv_to<arma::Mat<uint16_t>>::from(pixelMap);
        hdf5File.writeArray("/Images", imageName, uintMap);
    }

    // send the imagette to the client connected via tcp connection

    if (sendImagettesToClient)
    {
        Log.info("ClosedLoopDetectorWithMappedPSF: attempt to send the imagette to the client");

        sendImagetteToClient(&pixelMap, exposureNr);

    }


    if (numRowsSmearingMap != 0)
    {
        // Clear the string stream and compose the smearing map name

        myStream.str(string());      // insert empty string
        myStream.clear();            // clear eof bit

        myStream << "smearingMap" << setfill('0') << setw(6) << exposureNr;
        string smearingMapName = myStream.str();


        // Add the smearing map to the "SmearingMaps" group

        if (!includeQuantisation)
        {
            // Write the float array to HDF5

            //hdf5File.writeArray("/SmearingMaps", smearingMapName, smearingMap);
        }
        else
        {
            if ((smearingMap.min() < 0) || (smearingMap.max() >= (1 << 16)))
            {
                throw ConfigurationException("Detector: quantisation was applied but smearing map values are not in [0, 2^16[");
            }

            // Convert the float matrix to an unsigned uint16_t matrix

            arma::Mat<uint16_t> uintMap = arma::conv_to<arma::Mat<uint16_t>>::from(smearingMap);
            //hdf5File.writeArray("/SmearingMaps", smearingMapName, uintMap);
        }
        
    }

   // Clear the string stream and compose the bias map name

    myStream.str(string());      // insert empty string
    myStream.clear();            // clear eof bit

    myStream << "biasMap" << setfill('0') << setw(6) << exposureNr;
    string biasMapName = myStream.str();

    // Add the bias map to the "BiasMaps" group

    if (!includeQuantisation)
    {
        // Write the float array to HDF5

        //hdf5File.writeArray("/BiasMapsLeft", biasMapName, biasMapLeft);
        //hdf5File.writeArray("/BiasMapsRight", biasMapName, biasMapRight);
    }
    else
    {
        if ((biasMapLeft.min() < 0) || (biasMapLeft.max() >= (1 << 16)))
        {
            throw ConfigurationException("ClosedLoopDetectorWithMappedPSF: quantisation was applied but pixel values in the left bias map are not in [0, 2^16[");
        }

        if ((biasMapRight.min() < 0) || (biasMapRight.max() >= (1 << 16)))
        {
            throw ConfigurationException("ClosedLoopDetectorWithMappedPSF: quantisation was applied but pixel values in the right bias map are not in [0,2^16[");
        }

        // Convert the float matrix to an unsigned uint16_t matrix

        arma::Mat<uint16_t> uintMap = arma::conv_to<arma::Mat<uint16_t>>::from(biasMapLeft);
        //hdf5File.writeArray("/BiasMapsLeft", biasMapName, uintMap);

        uintMap = arma::conv_to<arma::Mat<uint16_t>>::from(biasMapRight);
        //hdf5File.writeArray("/BiasMapsRight", biasMapName, uintMap);
    }
    

    // Clear the string stream and compose the throughput map name

    myStream.str(string());      // insert empty string
    myStream.clear();            // clear eof bit

    myStream << "throughputMap" << setfill('0') << setw(6) << exposureNr;
    string throughputMapName = myStream.str();

    // Add the throughput map to the "ThroughputMaps" group

    //hdf5File.writeArray("/ThroughputMaps", throughputMapName, throughputMap);
}


 /**
 * \brief: variant of the generateFlatfieldMap without the writing of the maps to the hdf5 file - this is temporarity
 */

void ClosedLoopDetectorWithAnalyticNonGaussianPSF::generateFlatfieldMap()
{

    Log.info("ClosedLoopDetectorWithAnalyticNonGaussianPSF: generating flatfield map.");

    // Random number generation

    mt19937 flatfieldGenerator(flatfieldSeed);
    normal_distribution<double> flatfieldDistribution(0.0, 1.0);

    // Double the dimensions (this is necessary because of the behaviour of the Fourier transforms)

    int Nrows = 2 * numRowsPixelMap;
    int Ncolumns = 2 * numColumnsPixelMap;

    arma::cx_fmat evenMap = arma::cx_fmat(Nrows, Ncolumns);

    for(unsigned int row = 0; row < Nrows; row++)
    {
        for(unsigned int column = 0; column < Ncolumns; column++)
        {
            // Fourier space: generate white noise and include 1/f dependency

            evenMap(row, column) = flatfieldDistribution(flatfieldGenerator) / (pow(row, 2) + std::pow(column, 2) + 1);
        }
    }

    // Take the real part of the inverse Fourier transform

    evenMap = arma::ifft2(evenMap);
    arma::fmat realMap = arma::real(evenMap);

    // Cut out the appropriate part

    unsigned int numRowsFlatfield = Nrows / 2;
    unsigned int numColumnsFlatfield = Ncolumns / 2;
    
    flatfieldMap(arma::span::all, arma::span::all) = realMap(arma::span(0, numRowsFlatfield - 1), arma::span(0, numColumnsFlatfield - 1));
    flatfieldMap.reshape(numRowsFlatfield * numColumnsFlatfield, 1);

    // Normalisation
    //  - divide by mean and subtract 1.0 -> mean = 0.0
    //  - scale such that std.dev. = flatfield RMS and mean = 0.0
    //  - add 1.0

    flatfieldMap /= arma::mean(flatfieldMap.col(0));
    flatfieldMap -= 1;
    double scale = flatfieldNoiseRMS / arma::stddev(flatfieldMap.col(0));
    flatfieldMap *= scale;
    flatfieldMap += 1;

    flatfieldMap.reshape(numRowsFlatfield, numColumnsFlatfield);

    // Write the result to the HDF5 output file

    //Log.debug("ClosedLoopDetectorWithAnalyticNonGaussianPSF: writing PRNU to HDF5");

    //hdf5File.writeArray("/Flatfield", "PRNU", flatfieldMap);
}




/**
 * \brief change the window position and adapt all relevant pixel maps, if there was a new message from the server
 *        
 */
void ClosedLoopDetectorWithAnalyticNonGaussianPSF::setNewWindowPosition(std::tuple<bool, uint, uint, uint, uint, double> windowPositionTuple)
{
    // check whether a position change is necessary
    if (get<0>(windowPositionTuple))
    {
        // change the position of the subfield

        numRowsPixelMap         = get<1>(windowPositionTuple);
        numColumnsPixelMap      = get<2>(windowPositionTuple);

        subFieldZeroPointRow    = get<3>(windowPositionTuple);
        subFieldZeroPointColumn = get<4>(windowPositionTuple);
        
        orientationAngle      = deg2rad(get<5>(windowPositionTuple));
        
        Log.info("ClosedLoopDetectorWithMappedPSF: Changed numRowsPixelMap to: " + to_string(numRowsPixelMap));
        Log.info("ClosedLoopDetectorWithMappedPSF: Changed numColumnsPixelMap to: " + to_string(numColumnsPixelMap));

        Log.info("ClosedLoopDetectorWithMappedPSF: Changed subFieldZeroPointColumn to: " + to_string(subFieldZeroPointColumn));
        Log.info("ClosedLoopDetectorWithMappedPSF: Changed subFieldZeroPointRow to: " + to_string(subFieldZeroPointRow));

        Log.info("ClosedLoopDetectorWithMappedPSF: Changed orientationAngle to: " + to_string(orientationAngle));
            

        pixelMap.set_size(numRowsPixelMap, numColumnsPixelMap);

        biasMapLeft.set_size(numRowsBiasMap, numColumnsBiasMap);

        biasMapRight.set_size(numRowsBiasMap, numColumnsBiasMap);

        smearingMap.set_size(numRowsSmearingMap, numColumnsPixelMap);

        throughputMap.set_size(numRowsPixelMap, numColumnsPixelMap);

        throughputMap.ones(numRowsPixelMap, numColumnsPixelMap);

        // If we are going to apply open-shutter smearing, we have to know which pixels are within
        // the FOV (relevant only in case of mechanical vignetting).  When mechanical vignetting is
        // disabled, all pixels of the detector are inside the FOV.
            
        if(includeOpenShutterSmearing)
        {
            // Mechanical vignetting map:
            //  - no mechanical vignetting: all pixels of the sub-field inside FOV -> all values set to one
            //  - mechanical vignetting: set value of the pixels in the sub-field outside FOV to zero (others should be one) -> on creation of the throughput map
            
            mechanicalVignettingMask.ones(numRowsPixelMap, numColumnsPixelMap);

            // Number of exposed rows in each column:
            // - no mechanical vignetting: all exposed rows inside FOV (numRows - firstRowExposed)
            // - mechanical vignetting: count the exposed rows (i.e. from firstRowExposed) that are inside FOV
            
            numExposedRowsInFOV.zeros(numColumnsPixelMap);

            if(!includeMechanicalVignetting)
            {
                numExposedRowsInFOV.fill(numRows - firstRowExposed);
            }

            // Allocate memory for the different maps

            flatfieldMap.ones(numRowsPixelMap, numColumnsPixelMap);

            if(includeFlatfield)
            {
                // Generate the flatfield map

                generateFlatfieldMap();
            }
        }
    }
}



// ------------------------------------------------------------------------------------------------------------- //





/**
 * \brief override the takeExposure function to include the receiving of windowPosition from server 
 *        
 */
double ClosedLoopDetectorWithAnalyticGaussianPSF::takeExposure(int exposureNr, double startTime, double exposureTime)
{
    // Advance the internal clock until the given start time

    internalTime = startTime;

    // check if there are new updates for the window position

    if (getWindowPositionFromServer)
    {
        // receive window postion
        std::tuple<bool, uint, uint, uint, uint, double> windowPositionTuple = getNewWindowPosition(exposureTime);

        // set the window position according to the new message
        setNewWindowPosition(windowPositionTuple);
    }

    // Integration of point sources and background, taking into account jitter + drift.

    Log.info("ClosedLoopDetectorWithMappedPSF: Integrating light for exposure " + to_string(exposureNr) + " with exposure time = " + to_string(exposureTime));

    integrateLight(exposureNr, startTime, exposureTime);

    // Include noise effects like readout noise, photon noise, full well saturation, etc.
    // Note: readOut() needs the exposure time to compute the open shutter smearing.

    Log.info("ClosedLoopDetectorWithMappedPSF: Adding noise effects to exposure " + to_string(exposureNr));

    readOut(exposureTime);

    // Write the CCD subfield, the bias map, and the smearing map to the HDF5 file

    Log.debug("ClosedLoopDetectorWithMappedPSF: Writing PixelMap, smearing map, and bias map #" + to_string(exposureNr) + " to HDF5 file.");

    writePixelMapsToHDF5(exposureNr);

    // Advance the internal clock

    internalTime += exposureTime + readoutTimeBeforeNextExposure;

    return internalTime;
}


/**
 * \brief override the writePixelMapsToHDF5 function to include the sending of the imagettes to client 
 *        
 */
void ClosedLoopDetectorWithAnalyticGaussianPSF::writePixelMapsToHDF5(int exposureNr)
{
    // Compose the image name

    stringstream myStream;
    myStream << "image" << setfill('0') << setw(6) << exposureNr;
    string imageName = myStream.str();

    // Add the image to the "Images" group

    if (!includeQuantisation)
    {
        // Write the float array to HDF5

        hdf5File.writeArray("/Images", imageName, pixelMap);
    }
    else
    {
        // Write the pixel maps as 2-byte (16 bit) unsigned short integers.
        // As a safety check, first check that the extrema of the map are indeed
        // within the boundaries of such a data type.
       
        if((pixelMap.min() < 0) || (pixelMap.max() >= (1 << 16)))
        {
            throw ConfigurationException("ClosedLoopDetectorWithMappedPSF: quantisation was applied but pixel map values are not in [0, 2^16[");
        }

        // Convert the float matrix to an unsigned uint16_t matrix

        arma::Mat<uint16_t> uintMap = arma::conv_to<arma::Mat<uint16_t>>::from(pixelMap);
        hdf5File.writeArray("/Images", imageName, uintMap);
    }

    // send the imagette to the client connected via tcp connection

    if (sendImagettesToClient)
    {
        Log.info("ClosedLoopDetectorWithMappedPSF: attempt to send the imagette to the client");

        sendImagetteToClient(&pixelMap, exposureNr);

    }


    if (numRowsSmearingMap != 0)
    {
        // Clear the string stream and compose the smearing map name

        myStream.str(string());      // insert empty string
        myStream.clear();            // clear eof bit

        myStream << "smearingMap" << setfill('0') << setw(6) << exposureNr;
        string smearingMapName = myStream.str();


        // Add the smearing map to the "SmearingMaps" group

        if (!includeQuantisation)
        {
            // Write the float array to HDF5

            //hdf5File.writeArray("/SmearingMaps", smearingMapName, smearingMap);
        }
        else
        {
            if ((smearingMap.min() < 0) || (smearingMap.max() >= (1 << 16)))
            {
                throw ConfigurationException("Detector: quantisation was applied but smearing map values are not in [0, 2^16[");
            }

            // Convert the float matrix to an unsigned uint16_t matrix

            arma::Mat<uint16_t> uintMap = arma::conv_to<arma::Mat<uint16_t>>::from(smearingMap);
            //hdf5File.writeArray("/SmearingMaps", smearingMapName, uintMap);
        }
        
    }

   // Clear the string stream and compose the bias map name

    myStream.str(string());      // insert empty string
    myStream.clear();            // clear eof bit

    myStream << "biasMap" << setfill('0') << setw(6) << exposureNr;
    string biasMapName = myStream.str();

    // Add the bias map to the "BiasMaps" group

    if (!includeQuantisation)
    {
        // Write the float array to HDF5

        //hdf5File.writeArray("/BiasMapsLeft", biasMapName, biasMapLeft);
        //hdf5File.writeArray("/BiasMapsRight", biasMapName, biasMapRight);
    }
    else
    {
        if ((biasMapLeft.min() < 0) || (biasMapLeft.max() >= (1 << 16)))
        {
            throw ConfigurationException("ClosedLoopDetectorWithMappedPSF: quantisation was applied but pixel values in the left bias map are not in [0, 2^16[");
        }

        if ((biasMapRight.min() < 0) || (biasMapRight.max() >= (1 << 16)))
        {
            throw ConfigurationException("ClosedLoopDetectorWithMappedPSF: quantisation was applied but pixel values in the right bias map are not in [0,2^16[");
        }

        // Convert the float matrix to an unsigned uint16_t matrix

        arma::Mat<uint16_t> uintMap = arma::conv_to<arma::Mat<uint16_t>>::from(biasMapLeft);
        //hdf5File.writeArray("/BiasMapsLeft", biasMapName, uintMap);

        uintMap = arma::conv_to<arma::Mat<uint16_t>>::from(biasMapRight);
        //hdf5File.writeArray("/BiasMapsRight", biasMapName, uintMap);
    }
    

    // Clear the string stream and compose the throughput map name

    myStream.str(string());      // insert empty string
    myStream.clear();            // clear eof bit

    myStream << "throughputMap" << setfill('0') << setw(6) << exposureNr;
    string throughputMapName = myStream.str();

    // Add the throughput map to the "ThroughputMaps" group

    //hdf5File.writeArray("/ThroughputMaps", throughputMapName, throughputMap);
}

 /**
 * \brief: variant of the generateFlatfieldMap without the writing of the maps to the hdf5 file - this is temporarity
 */
void ClosedLoopDetectorWithAnalyticGaussianPSF::generateFlatfieldMap()
{

    Log.info("ClosedLoopDetectorWithAnalyticGaussianPSF: generating flatfield map.");

    // Random number generation

    mt19937 flatfieldGenerator(flatfieldSeed);
    normal_distribution<double> flatfieldDistribution(0.0, 1.0);

    // Double the dimensions (this is necessary because of the behaviour of the Fourier transforms)

    int Nrows = 2 * numRowsPixelMap;
    int Ncolumns = 2 * numColumnsPixelMap;

    arma::cx_fmat evenMap = arma::cx_fmat(Nrows, Ncolumns);

    for(unsigned int row = 0; row < Nrows; row++)
    {
        for(unsigned int column = 0; column < Ncolumns; column++)
        {
            // Fourier space: generate white noise and include 1/f dependency
            // (Note: see https://en.wikipedia.org/wiki/Pink_noise#Generalization_to_more_than_one_dimension)

            evenMap(row, column) = flatfieldDistribution(flatfieldGenerator) / (pow(row, 2) + std::pow(column, 2) + 1);
        }
    }

    // Take the real part of the inverse Fourier transform

    evenMap = arma::ifft2(evenMap);
    arma::fmat realMap = arma::real(evenMap);

    // Cut out the appropriate part

    unsigned int numRowsFlatfield = Nrows / 2;
    unsigned int numColumnsFlatfield = Ncolumns / 2;
    
    flatfieldMap(arma::span::all, arma::span::all) = realMap(arma::span(0, numRowsFlatfield - 1), arma::span(0, numColumnsFlatfield - 1));
    flatfieldMap.reshape(numRowsFlatfield * numColumnsFlatfield, 1);

    // Normalisation
    //  - divide by mean and subtract 1.0 -> mean = 0.0
    //  - scale such that std.dev. = flatfield RMS and mean = 0.0
    //  - add 1.0

    flatfieldMap /= arma::mean(flatfieldMap.col(0));
    flatfieldMap -= 1;
    double scale = flatfieldNoiseRMS / arma::stddev(flatfieldMap.col(0));
    flatfieldMap *= scale;
    flatfieldMap += 1;

    flatfieldMap.reshape(numRowsFlatfield, numColumnsFlatfield);

    // Write the result to the HDF5 output file

/*    Log.debug("Detector: writing PRNU to HDF5");

    hdf5File.writeArray("/Flatfield", "PRNU", flatfieldMap);*/
}






/**
 * \brief change the window position and adapt all relevant pixel maps, if there was a new message from the server
 *        
 */
void ClosedLoopDetectorWithAnalyticGaussianPSF::setNewWindowPosition(std::tuple<bool, uint, uint, uint, uint, double> windowPositionTuple)
{
    // check whether a position change is necessary
    if (get<0>(windowPositionTuple))
    {
        // change the position of the subfield

        numRowsPixelMap         = get<1>(windowPositionTuple);
        numColumnsPixelMap      = get<2>(windowPositionTuple);

        subFieldZeroPointRow    = get<3>(windowPositionTuple);
        subFieldZeroPointColumn = get<4>(windowPositionTuple);
        
        orientationAngle      = deg2rad(get<5>(windowPositionTuple));
        
        Log.info("ClosedLoopDetectorWithMappedPSF: Changed numRowsPixelMap to: " + to_string(numRowsPixelMap));
        Log.info("ClosedLoopDetectorWithMappedPSF: Changed numColumnsPixelMap to: " + to_string(numColumnsPixelMap));

        Log.info("ClosedLoopDetectorWithMappedPSF: Changed subFieldZeroPointColumn to: " + to_string(subFieldZeroPointColumn));
        Log.info("ClosedLoopDetectorWithMappedPSF: Changed subFieldZeroPointRow to: " + to_string(subFieldZeroPointRow));

        Log.info("ClosedLoopDetectorWithMappedPSF: Changed orientationAngle to: " + to_string(orientationAngle));
            

        pixelMap.set_size(numRowsPixelMap, numColumnsPixelMap);

        biasMapLeft.set_size(numRowsBiasMap, numColumnsBiasMap);

        biasMapRight.set_size(numRowsBiasMap, numColumnsBiasMap);

        smearingMap.set_size(numRowsSmearingMap, numColumnsPixelMap);

        throughputMap.set_size(numRowsPixelMap, numColumnsPixelMap);

        throughputMap.ones(numRowsPixelMap, numColumnsPixelMap);

        // If we are going to apply open-shutter smearing, we have to know which pixels are within
        // the FOV (relevant only in case of mechanical vignetting).  When mechanical vignetting is
        // disabled, all pixels of the detector are inside the FOV.
            
        if(includeOpenShutterSmearing)
        {
            // Mechanical vignetting map:
            //  - no mechanical vignetting: all pixels of the sub-field inside FOV -> all values set to one
            //  - mechanical vignetting: set value of the pixels in the sub-field outside FOV to zero (others should be one) -> on creation of the throughput map
            
            mechanicalVignettingMask.ones(numRowsPixelMap, numColumnsPixelMap);

            // Number of exposed rows in each column:
            // - no mechanical vignetting: all exposed rows inside FOV (numRows - firstRowExposed)
            // - mechanical vignetting: count the exposed rows (i.e. from firstRowExposed) that are inside FOV
            
            numExposedRowsInFOV.zeros(numColumnsPixelMap);

            if(!includeMechanicalVignetting)
            {
                numExposedRowsInFOV.fill(numRows - firstRowExposed);
            }

            // Allocate memory for the different maps

            flatfieldMap.ones(numRowsPixelMap, numColumnsPixelMap);

            if(includeFlatfield)
            {
                // Generate the flatfield map

                generateFlatfieldMap();
            }
        }
    }
}

