#include "Photometry.h"




/**
 * \brief      Constructor
 *
 * \param configParam   Configuration parameters for the Camera
 * \param hdf5file      HFD5 file to write the camera information to
 * \param telescope     Telescope on which the camera is mounted on
 * \param sky           Sky object to query which stars we are seeing
 */

Photometry::Photometry(ConfigurationParameters &configParam,
		       HDF5File &hdf5file,
		       Camera &camera)
  : HDF5Writer(hdf5file)
{
    // Parse the parameters from the configuration file.

    configure(configParam);
}




    
/**
 * \brief Preprocessing of the subfield before photometric extraction can take place.
 *
 */

void Photometry::Preprocessing(const unsigned int exposureNr)
{
    const unsigned int zeroBasedExposureNr = exposureNr - beginExposureNr;
    const double varianceRON = sqrt(pow(readoutNoise, 2) + pow(frontEndElectronics->getReadoutNoise(), 2));      // [electrons / pixel]
 
    // Make a (deep) copy of the pixelMap on which we can do some reductions without altering the original pixelMap

    arma::Mat<float> image(pixelMap);

    // Subtract the bias

    float meanBias = 0.0;
    if (subFieldZeroPointColumn <  numColumns / 2)
    {
        meanBias = arma::mean(arma::mean(biasMapLeft));
    }
    else
    {
        meanBias = arma::mean(arma::mean(biasMapRight));
    }
    image -= meanBias;

    // Correct for open-shutter smearing

    image.each_row() -= arma::mean(smearingMap - meanBias, 0);

    // Convert from [ADU] to [electrons] using the gain

    if (subFieldZeroPointColumn <  numColumns / 2)
    {
        image /= combinedGainLeft;
    }
    else
    {
        image /= combinedGainRight;
    }
    
    // Subtract the sky background

    const double skyBackground = camera.getTotalSkyBackground();  // [photons/pixel/exposure]
    image -= throughputMap * skyBackground;                       // [e-/pixel/exposure]
}







/**
 * \brief Extract the photometric light curve for a specified list of stars
 *
 * TODO: - better error catching when the stars for which a lightcurve is requested are (sometimes) not in the subfield
 *       - better treatment when there are no contaminants
 */

void Photometry::extractPhotometry(const unsigned int exposureNr)
{

    // Nr of stars for which we want a lightcurve
  
    const int Ntargets = photStarIDs.size();
    if (Ntargets == 0)
    {
        Log.warning("Photometry:extractPhotometry: no stars found for which photometry is requested. Skipping extractPhotometry()."); 
        return;
    }

    // Loop over all targets for which you need a lightcurve
    
    for (int n = 0; n < Ntargets; n++)
    {
        // Collect info on the position and the input flux of the target

        int starID = photStarIDs[n];

        double time;         // Time stamp of the last exposure         [s]
        double xFPtarget;    // Mean x-coordinate in the focal plane    [mm]
        double yFPtarget;    // Mean y-coordinate in the focal plane    [mm] 
        double rowTarget;    // Mean row coordinate in the subfield     [pix] 
        double colTarget;    // Mean column coordinate in the subfield  [pix]
        double fluxTarget;   // Total flux during the exposure          [photons/exposure]
    
        tie(time, xFPtarget, yFPtarget, rowTarget, colTarget, fluxTarget) = camera.getInfoForTheMostRecentExposureForStar(starID);

        if (fluxTarget == -1.0)
        {
            Log.warning("Photometry:extractPhotometry: no info found for star " + to_string(starID) + " for which photometry is requested");
            continue;
        }

        inputFluxTarget.at(starID).at(zeroBasedExposureNr) = fluxTarget;

        // If this is the first exposure, the mask is alway defined for the current target

        double timeSinceLastMaskUpdate = 0.0;

        if (exposureNr != beginExposureNr)
        {
            timeSinceLastMaskUpdate = (exposureNr - exposureNrOfMaskUpdate.at(starID).back()) * cycleTime;  // [s]
        }

        if ((exposureNr == beginExposureNr) || (timeSinceLastMaskUpdate > maskUpdateInterval))
        {
            Log.debug("Photometry::extractPhotometry: updating mask of star ID " + to_string(starID) + " for exposure " + to_string(exposureNr));
            Log.debug("Photometry::extractPhotometry: creating single-target and contamination maps");

            arma::Mat<float> singleTargetMap(numRowsPixelMap, numColumnsPixelMap);
            arma::Mat<float> contaminantMap(numRowsPixelMap, numColumnsPixelMap);

            // Create a noiseless subfield as if there was only the flux of this single target

            singleTargetMap.zeros();
            double r = rad2deg(camera.getGnomonicRadialDistanceFromOpticalAxis(xFPtarget, yFPtarget));
            double p = atan2(yFPtarget, xFPtarget);
            bool success = addFluxToMap(singleTargetMap, rowTarget, colTarget, r, p, fluxTarget);

            // Create a noiseless subfield of only the possible contaminants

            contaminantMap.zeros();

            int numContaminants = 0;
            starInfoIterator begin, end;
            tie(begin, end) = camera.getInfoForTheMostRecentExposureForAllStars();
            for (auto it = begin; it != end; it++)
            {
                if (it->first == starID) continue;                        // A star is never its own contaminant 
                double xFPcont =  (it->second)[0] / (it->second)[5];      // [mm]
                double yFPcont =  (it->second)[1] / (it->second)[5];      // [mm]
                double rowCont =  (it->second)[2] / (it->second)[5];      // [pix]
                double colCont =  (it->second)[3] / (it->second)[5];      // [pix]
                double fluxCont = (it->second)[4];                        // [photons/exposure]

                // Skip the contaminants that are too distant from the target to have any effect

                if ((abs(colCont - colTarget) > contaminationRadius) or (abs(rowCont - rowTarget) > contaminationRadius))
                    continue;
                
                // Add the PSF of the contaminant to the contaminant map
                
                r = rad2deg(camera.getGnomonicRadialDistanceFromOpticalAxis(xFPcont, yFPcont));
                p = atan2(yFPcont, xFPcont);
                success = addFluxToMap(contaminantMap, rowCont, colCont, r, p, fluxCont);
                numContaminants++;
            }

            Log.debug("Photometry::extractPhotometry: Found " + to_string(numContaminants) + " contaminants for star ID " + to_string(starID));
            Log.debug("Photometry::extractPhotometry: selecting which pixels belong to the mask");

            // For the mask of our target will only consider a 7x7 area around the barycenter. Get the boundaries of that area.

            const int minRow = max(0, int(rowTarget)-3);
            const int maxRow = min(int(numRowsPixelMap) - 1, int(rowTarget)+3);                         // maxRow inclusive
            const int minCol = max(0, int(colTarget)-3);
            const int maxCol = min(int(numColumnsPixelMap) - 1, int(colTarget)+3);                      // maxCol inclusive
            
            Log.debug("Photometry::extractPhotometry: determining mask within the area: pixelMap rows: "
                      + to_string(minRow) + " -> " + to_string(maxRow) + ", cols: "
                      + to_string(minCol) + " -> " + to_string(maxCol) + ". End points inclusive");
           
            if ((numRowsPixelMap <= 7) || (numColumnsPixelMap <= 7))
            {
                Log.warning("Photometry::extractPhotometry: size of pixel map is smaller than 8x8");
            }

            // For the pixels in the designated area around our target, compute the variance and the noise/signal ratio of the signal.
            // Example size: if the pixelMap is 100x100 pixels, and we consider a mask of 4x4 pixels, then NSRmap is a 2D array of size
            //               100x100, but flatNSRmap is a 1D array of size 16. 

            arma::Mat<float> NSRmap(numRowsPixelMap, numColumnsPixelMap, arma::fill::zeros); 
            arma::Mat<float> varianceMap(numRowsPixelMap, numColumnsPixelMap, arma::fill::zeros);

            vector<double> flatNSRmap;   
            for (int irow = minRow; irow <= maxRow; irow++)
            {
                for (int icol = minCol; icol <= maxCol; icol++)
                {
                    // We assume photon noise, so the variance equals the flux. We multiply by the throughput so that both terms
                    // are expressed in [e-/exposure]. 

                    varianceMap(irow, icol) = (singleTargetMap(irow, icol) + contaminantMap(irow, icol) + skyBackground) * throughputMap(irow, icol) + varianceRON;
                    NSRmap(irow, icol) = sqrt(varianceMap(irow, icol)) / singleTargetMap(irow, icol); 
                    flatNSRmap.push_back(NSRmap(irow, icol));
                }
            }

            // Order the pixels in the (flattened) NSR map from low to high N/S ratio (i.e. high to low S/N)

            vector<unsigned int> indices(flatNSRmap.size());
            iota(indices.begin(), indices.end(), 0);
            stable_sort(indices.begin(), indices.end(), [&flatNSRmap](unsigned int i, unsigned int j) {return flatNSRmap[i] < flatNSRmap[j];});
            vector<unsigned int> rowIndex(flatNSRmap.size());           
            vector<unsigned int> colIndex(flatNSRmap.size());
            const int NcolsMask = maxCol-minCol +1; 

	    // Transform from indices in flatNSRmap to indices in NSRmap
	    
            for (int i = 0; i < rowIndex.size(); i++)
            {
                rowIndex[i] = minRow + (unsigned int)(indices[i]) / NcolsMask;
                colIndex[i] = minCol + (unsigned int)(indices[i]) % NcolsMask; 
            }

            // Build the mask, starting with the pixel with the lowest NSR, adding one pixel at the time,
            // with the condition that adding a pixel should contribute more to the aggregated signal than
	    // to the aggregated noise. Initialize with the first pixel

            double aggregatedVariance            = varianceMap(rowIndex[0], colIndex[0]);
            double aggregatedSingleTargetFlux    = singleTargetMap(rowIndex[0], colIndex[0]);
            double aggregatedObservedTargetFlux  = image(rowIndex[0], colIndex[0]);
            double aggregatedNSR                 = NSRmap(rowIndex[0], colIndex[0]);
            maskSizeTarget[starID].push_back(1);
 
            rowIndexOfMaskOfTarget[starID][exposureNr] = {rowIndex[0]}; 
            colIndexOfMaskOfTarget[starID][exposureNr] = {colIndex[0]};

            // Then add other pixels

            for (int i = 1; i < rowIndex.size(); i++)
            {
                double temp = sqrt(aggregatedVariance + varianceMap(rowIndex[i], colIndex[i])) / (aggregatedSingleTargetFlux + singleTargetMap(rowIndex[i], colIndex[i]));
                if (temp < aggregatedNSR)
                {
                    // The aggregated NSR improved by adding a pixel, so include the pixel in the mask

                    aggregatedVariance           += varianceMap(rowIndex[i], colIndex[i]);
                    aggregatedSingleTargetFlux   += singleTargetMap(rowIndex[i], colIndex[i]);
                    aggregatedObservedTargetFlux += image(rowIndex[i], colIndex[i]);
                    aggregatedNSR = temp;
                    maskSizeTarget.at(starID).back() += 1;
                    rowIndexOfMaskOfTarget.at(starID).at(exposureNr).push_back(rowIndex[i]);
                    colIndexOfMaskOfTarget.at(starID).at(exposureNr).push_back(colIndex[i]);
                }
                else
                {
                    // The aggregated NSR did not improve by adding this pixel. Not only can we ignore exclude this pixel from the 
                    // mask, but also all subsequent ones that have an even worse noise/signal ratio. So finalize the mask for this target, and 
                    // then break out of the for-loop.

                    estimatedFluxTarget.at(starID).at(zeroBasedExposureNr) = aggregatedObservedTargetFlux; 
                    varFluxTarget.at(starID).at(zeroBasedExposureNr) = aggregatedVariance; 
                    NSRtarget.at(starID).push_back(aggregatedNSR);

                    // Disregard all other pixels of the window around the target star: they all contribute more to the noise than to the signal.

                    break;
                }
            }

            // Update the exposure nr of this mask update to the current exposure number

            exposureNrOfMaskUpdate.at(starID).push_back(exposureNr);

        }
        else
        {
            // For all other exposures, simply use the same (most recent) mask. We reuse the NSR of the previous mask, so no need to recompute it.

            Log.debug("Detector::applyPhotometry: extracting flux for target ID " + to_string(starID) + " for exposure " + to_string(exposureNr) + " with an old mask");

            const unsigned int exposureNrOfLastMaskUpdate =  exposureNrOfMaskUpdate.at(starID).back();
            estimatedFluxTarget.at(starID).at(zeroBasedExposureNr) = 0.0;
            for (int j = 0; j < maskSizeTarget[starID].back(); j++)
            {
                const unsigned int rowIndex = rowIndexOfMaskOfTarget.at(starID).at(exposureNrOfLastMaskUpdate).at(j);
                const unsigned int colIndex = colIndexOfMaskOfTarget.at(starID).at(exposureNrOfLastMaskUpdate).at(j);
                estimatedFluxTarget.at(starID).at(zeroBasedExposureNr) += image(rowIndex, colIndex);
            }
        }
    } // end loop over all targets for which we want light curves
} // end applyPhotometry()








/**
 * \brief Apply photometry: first preprocessing and then extract photometry.
 *
 */

void Photometry::applyPhotometry(const unsigned int exposureNr)
{

  // Preprocessing

  preprocessing(exposureNr);

  // Extract photometry

  extractPhotometry(exposureNr);
    
}








/**
 * \brief Apply photometry
 *
 */

void Photometry::writePhotometry()
{

    Log.info("Writing photometry to the HDF5 file");

    // Create major HDF5 directories
    
    hdf5File.createGroup("/Photometry");
    hdf5File.createGroup("/Photometry/Masks");
    hdf5File.createGroup("/Photometry/Lightcurves");

    // Since mask updates are done for all stars we pick the first star ID
    
    string groupName = "/Photometry/Masks";
    string arrayName = "exposureNrOfMaskUpdate";
    int starID = photStarIDs[0];
    hdf5File.writeArray(groupName, arrayName, exposureNrOfMaskUpdate[starID].data(), exposureNrOfMaskUpdate[starID].size());

    for (auto starID : photStarIDs)
    {
	string starName = to_string(starID);
	groupName = "/Photometry/Lightcurves/starID" + starName;
	hdf5File.createGroup(groupName);

	arrayName = "inputFlux";
	hdf5File.writeArray(groupName, arrayName, inputFluxTarget[starID].data(), inputFluxTarget[starID].size());

	arrayName = "estimatedFlux";
	hdf5File.writeArray(groupName, arrayName, estimatedFluxTarget[starID].data(), estimatedFluxTarget[starID].size());

	groupName = "/Photometry/Masks/starID" + starName;
	hdf5File.createGroup(groupName);

	arrayName = "maskSize";
	hdf5File.writeArray(groupName, arrayName, maskSizeTarget[starID].data(), maskSizeTarget[starID].size());

	arrayName = "maskNSR";
	hdf5File.writeArray(groupName, arrayName, NSRtarget[starID].data(), NSRtarget[starID].size());

	for(auto iter = rowIndexOfMaskOfTarget[starID].begin(); iter != rowIndexOfMaskOfTarget[starID].end(); ++iter)
	{
	    const unsigned int exposureNumber = iter->first;

	    stringstream myStream;
	    myStream << "Exposure" << setfill('0') << setw(6) << exposureNumber;
	    groupName = "/Photometry/Masks/starID" + starName + "/" + myStream.str();
	    hdf5File.createGroup(groupName);

	    arrayName = "maskRowIndices"; 
	    hdf5File.writeArray(groupName, arrayName, rowIndexOfMaskOfTarget[starID][exposureNumber].data(), rowIndexOfMaskOfTarget[starID][exposureNumber].size());

	    arrayName = "maskColumnIndices";
	    hdf5File.writeArray(groupName, arrayName, colIndexOfMaskOfTarget[starID][exposureNumber].data(), colIndexOfMaskOfTarget[starID][exposureNumber].size());
            }
        }
}
