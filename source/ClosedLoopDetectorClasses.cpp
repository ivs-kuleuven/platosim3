#include "ClosedLoopDetectorClasses.h"


/**
 * \brief override the takeExposure function to include the receiving of windowPosition from server 
 *        
 */
double ClosedLoopDetectorWithSymmetricalMappedPSF::takeExposure(int exposureNr, double startTime, double exposureTime)
{

    // check if there are new updates for the window position

    if (getWindowPositionFromServer)
    {
        // receive window postion
        std::tuple<bool, uint, uint, uint, uint, double> windowPositionTuple = getNewWindowPosition(exposureTime);

        // set the window position according to the new message
        setNewWindowPosition(windowPositionTuple);
    }


    DetectorWithSymmetricalMappedPSF::takeExposure(exposureNr, startTime, exposureTime);
}


/**
 * \brief override the writePixelMapsToHDF5 function to include the sending of the imagettes to client 
 *        
 */
void ClosedLoopDetectorWithSymmetricalMappedPSF::writePixelMapsToHDF5(int exposureNr)
{
    DetectorWithSymmetricalMappedPSF::writePixelMapsToHDF5(exposureNr);


    if (sendImagettesToClient)
    {
        Log.info("ClosedLoopDetectorWithSymmetricalMappedPSF: attempt to send the imagette to the client");

        sendImagetteToClient(&pixelMap, exposureNr);

    }

}


/**
 * \brief change the window position and adapt all relevant pixel maps, if there was a new message from the server
 *        
 */
void ClosedLoopDetectorWithSymmetricalMappedPSF::setNewWindowPosition(std::tuple<bool, uint, uint, uint, uint, double> windowPositionTuple)
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
        
        Log.info("ClosedLoopDetectorWithSymmetricalMappedPSF: Changed numRowsPixelMap to: " + to_string(numRowsPixelMap));
        Log.info("ClosedLoopDetectorWithSymmetricalMappedPSF: Changed numColumnsPixelMap to: " + to_string(numColumnsPixelMap));

        Log.info("ClosedLoopDetectorWithSymmetricalMappedPSF: Changed subFieldZeroPointColumn to: " + to_string(subFieldZeroPointColumn));
        Log.info("ClosedLoopDetectorWithSymmetricalMappedPSF: Changed subFieldZeroPointRow to: " + to_string(subFieldZeroPointRow));

        Log.info("ClosedLoopDetectorWithSymmetricalMappedPSF: Changed orientationAngle to: " + to_string(orientationAngle));
            

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
                for (int subsubfieldx = 0; subsubfieldx < numsubsubfieldsx; subsubfieldx++)  //%% Loop to select a different psf for each subsubfield
    {
        for (int subsubfieldy = 0; subsubfieldy < numsubsubfieldsy; subsubfieldy++)
        {
    setPsfForSubfield(subsubfieldx, subsubfieldy);  //%% Added subsubfield
        }
    }

        }
    }
}



// ------------------------------------------------------------------------------------------------------------- //




/**
 * \brief override the takeExposure function to include the receiving of windowPosition from server 
 *        
 */
double ClosedLoopDetectorWithAsymmetricalMappedPSF::takeExposure(int exposureNr, double startTime, double exposureTime)
{

    // check if there are new updates for the window position

    if (getWindowPositionFromServer)
    {
        // receive window postion
        std::tuple<bool, uint, uint, uint, uint, double> windowPositionTuple = getNewWindowPosition(exposureTime);

        // set the window position according to the new message
        setNewWindowPosition(windowPositionTuple);
    }


    DetectorWithAsymmetricalMappedPSF::takeExposure(exposureNr, startTime, exposureTime);
}


/**
 * \brief override the writePixelMapsToHDF5 function to include the sending of the imagettes to client 
 *        
 */
void ClosedLoopDetectorWithAsymmetricalMappedPSF::writePixelMapsToHDF5(int exposureNr)
{
    DetectorWithAsymmetricalMappedPSF::writePixelMapsToHDF5(exposureNr);


    if (sendImagettesToClient)
    {
        Log.info("ClosedLoopDetectorWithAsymmetricalMappedPSF: attempt to send the imagette to the client");

        sendImagetteToClient(&pixelMap, exposureNr);

    }

}


/**
 * \brief change the window position and adapt all relevant pixel maps, if there was a new message from the server
 *        
 */
void ClosedLoopDetectorWithAsymmetricalMappedPSF::setNewWindowPosition(std::tuple<bool, uint, uint, uint, uint, double> windowPositionTuple)
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
        
        Log.info("ClosedLoopDetectorWithAsymmetricalMappedPSF: Changed numRowsPixelMap to: " + to_string(numRowsPixelMap));
        Log.info("ClosedLoopDetectorWithAsymmetricalMappedPSF: Changed numColumnsPixelMap to: " + to_string(numColumnsPixelMap));

        Log.info("ClosedLoopDetectorWithAsymmetricalMappedPSF: Changed subFieldZeroPointColumn to: " + to_string(subFieldZeroPointColumn));
        Log.info("ClosedLoopDetectorWithAsymmetricalMappedPSF: Changed subFieldZeroPointRow to: " + to_string(subFieldZeroPointRow));

        Log.info("ClosedLoopDetectorWithAsymmetricalMappedPSF: Changed orientationAngle to: " + to_string(orientationAngle));
            

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
             for (int subsubfieldx = 0; subsubfieldx < numsubsubfieldsx; subsubfieldx++)  //%% Loop to select a different psf for each subsubfield
    {
        for (int subsubfieldy = 0; subsubfieldy < numsubsubfieldsy; subsubfieldy++)
        {
    setPsfForSubfield(subsubfieldx, subsubfieldy);  //%% Added subsubfield
        }
    }

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

    // check if there are new updates for the window position

    if (getWindowPositionFromServer)
    {
        // receive window postion
        std::tuple<bool, uint, uint, uint, uint, double> windowPositionTuple = getNewWindowPosition(exposureTime);

        // set the window position according to the new message
        setNewWindowPosition(windowPositionTuple);
    }


    DetectorWithAnalyticNonGaussianPSF::takeExposure(exposureNr, startTime, exposureTime);
}


/**
 * \brief override the writePixelMapsToHDF5 function to include the sending of the imagettes to client 
 *        
 */
void ClosedLoopDetectorWithAnalyticNonGaussianPSF::writePixelMapsToHDF5(int exposureNr)
{
    DetectorWithAnalyticNonGaussianPSF::writePixelMapsToHDF5(exposureNr);


    if (sendImagettesToClient)
    {
        Log.info("ClosedLoopDetectorWithAnalyticNonGaussianPSF: attempt to send the imagette to the client");

        sendImagetteToClient(&pixelMap, exposureNr);

    }
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

    // check if there are new updates for the window position

    if (getWindowPositionFromServer)
    {
        // receive window postion
        std::tuple<bool, uint, uint, uint, uint, double> windowPositionTuple = getNewWindowPosition(exposureTime);

        // set the window position according to the new message
        setNewWindowPosition(windowPositionTuple);
    }

    DetectorWithAnalyticGaussianPSF::takeExposure(exposureNr, startTime, exposureTime);
}


/**
 * \brief override the writePixelMapsToHDF5 function to include the sending of the imagettes to client 
 *        
 */
void ClosedLoopDetectorWithAnalyticGaussianPSF::writePixelMapsToHDF5(int exposureNr)
{
    DetectorWithAnalyticGaussianPSF::writePixelMapsToHDF5(exposureNr);


    if (sendImagettesToClient)
    {
        Log.info("ClosedLoopDetectorWithAnalyticGaussianPSF: attempt to send the imagette to the client");

        sendImagetteToClient(&pixelMap, exposureNr);

    }
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

