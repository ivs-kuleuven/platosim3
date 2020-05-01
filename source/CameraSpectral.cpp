#include "CameraSpectral.h"

/**
 * \brief      Constructor
 *
 * \param configParam   Configuration parameters for the Camera
 * \param hdf5file      HFD5 file to write the camera information to
 * \param telescope     Telescope on which the camera is mounted on
 * \param sky           Sky object to query which stars we are seeing
 */

CameraSpectral::CameraSpectral(ConfigurationParameters &configParam, HDF5File &hdf5file, Platform &platform, Telescope &telescope, Sky &sky) : Camera(configParam, hdf5file, platform, telescope, sky)
{
    Spectral = new SpectralDependenceUtility(configParam);
    getParameters(configParam);
}



/**
 *
 */
void CameraSpectral::getParameters(ConfigurationParameters &configParam)
{
    //SpectralDependenceUtility Spectral(configParam);
    binnumber = Spectral->binnumber;
    binwidth = Spectral->binwidth;
    lowerWavelength = Spectral->lowerWavelength;
    referenceWavelength = Spectral->referenceWavelength;
    useQE = Spectral->useQE;
    if (useQE)
    {
        QESpectral = Spectral->QESpectral;
        meanQE = Spectral->meanQE;
    }
    else
    {
        QESpectral.assign(binnumber,1.); 
        meanQE = 1.;
    }
}


/**
 *
 */
void CameraSpectral::exposeDetector(Detector &detector, double startTime, double exposureTime, double readoutTimeBeforeNextExposure)
{
    // Get the value for the degrading TransmissionEfficiency parameter at the startTime of this exposure

    double transmissionEfficiency = telescope.getTransmissionEfficiency(startTime);

    Log.debug("Camera: TransmissionEfficiency at time "+to_string(startTime)+" is "+to_string(transmissionEfficiency));

    vector<double> transmissionEfficiencySpectral = Spectral->getSpectralTransmissionEfficiency(startTime);


    // Get the focal plane coordinates of the center and the corners of the subfield (in [mm]).
    // To compute the diagonal length of the subfield, we only need the lower left (X00, Y00)
    // and the upper right (X11, Y11) corner of the subfield.

    double centerXmm, centerYmm;
    tie(centerXmm, centerYmm) = detector.getFocalPlaneCoordinatesOfSubfieldCenter();

    double corner00Xmm, corner00Ymm, corner11Xmm, corner11Ymm, dummy;
    tie(corner00Xmm, corner00Ymm, dummy, dummy, corner11Xmm, corner11Ymm, dummy, dummy) = detector.getFocalPlaneCoordinatesOfSubfieldCorners();

    // Convert the undistorted [mm] to distorted [mm] focal plane coordinates

    if (includeFieldDistortion)
    {
        Log.info("Camera: including field distortion");

        tie(centerXmm, centerYmm) = distortedToUndistortedFocalPlaneCoordinates(centerXmm, centerYmm);
        tie(corner00Xmm, corner00Ymm) = distortedToUndistortedFocalPlaneCoordinates(corner00Xmm, corner00Ymm);
        tie(corner11Xmm, corner11Ymm) = distortedToUndistortedFocalPlaneCoordinates(corner11Xmm, corner11Ymm);
    }

    Log.debug("Camera: center of subfield at (Xmm, Ymm) = (" + to_string(centerXmm) + ", " + to_string(centerYmm) + ") mm");
    Log.debug("Camera: lower left corner of subfield at (Xmm, Ymm) = (" + to_string(corner00Xmm) + ", " + to_string(corner00Ymm) + ") mm");
    Log.debug("Camera: upper right corner of subfield at (Xmm, Ymm) = (" + to_string(corner11Xmm) + ", " + to_string(corner11Ymm) + ") mm");


    // Convert the focal plane coordinates [mm] to (alpha, delta) equatorial sky coordinates [rad]

    double centerRA, centerDec;
    tie(centerRA, centerDec) = focalPlaneToSkyCoordinates(centerXmm, centerYmm);

    double corner00RA, corner00Dec;
    tie(corner00RA, corner00Dec) = focalPlaneToSkyCoordinates(corner00Xmm, corner00Ymm);

    double corner11RA, corner11Dec;
    tie(corner11RA, corner11Dec) = focalPlaneToSkyCoordinates(corner11Xmm, corner11Ymm);

    Log.debug("Camera: center of subfield at (alpha, delta) = (" + to_string(rad2deg(centerRA)) + ", " + to_string(rad2deg(centerDec)) + ") deg");
    Log.debug("Camera: lower left corner of subfield at (alpha, delta) = (" + to_string(rad2deg(corner00RA)) + ", " + to_string(rad2deg(corner00Dec)) + ") deg");
    Log.debug("Camera: upper right corner of subfield at (alpha, delta) = (" + to_string(rad2deg(corner11RA)) + ", " + to_string(rad2deg(corner11Dec)) + ") deg");


    // Compute the angular distance on the sky between the lower left and the upper right corner
    // of the subfield, to estimate the "radius" of the subfield.

    SkyCoordinates skyCoordinates00(corner00RA, corner00Dec, Angle::radians);
    SkyCoordinates skyCoordinates11(corner11RA, corner11Dec, Angle::radians);
    
    double radius = angularDistanceBetween(skyCoordinates00, skyCoordinates11, Angle::radians) / 2.0;

    Log.debug("Camera: semi-diagonal of subfield = " + to_string(rad2deg(radius)) + " deg");


    // Get a catalog of stars that fall on the subfield. Take the radius a bit larger so that the 
    // queried area includes possible small shifts of the projected subfield because of jitter.

    const unsigned long Nstars = sky.selectStarsWithinRadiusFrom(centerRA, centerDec, radius * 1.1, Angle::radians);

    Log.info("Camera: Found " + to_string(Nstars) + " stars on and near the subfield");  

    if (includeAberrationCorrection)
    {
        Log.info("Camera: applying " + aberrationCorrectionType + " aberration correction to the selected stars in the subfield.");

        // The time at the middle of the time series is the time when the Sun is defined to be 180 degrees away from platform pointing

        double timeMiddle = numExposures * (exposureTime + readoutTimeBeforeNextExposure) / 2.0;

        // Get the apparent position of the stars, i.e. apply the differential aberration correction to
        // all the star positions in this starCatalog.

        // We do this calcuation only once per exposure as the effect is negligible within the exposure time
    
        sky.aberrateSelectedStarPositions(platform, aberrationCorrectionType, startTime, timeMiddle);
    }


    // If the telescope and/or platform show small variations (e.g. due to jitter) during the exposure,
    // the exposure time is split up in many small intervals, to track the effect of these variations
    // on the exposure. The largest time interval for which the variations can still be reliably tracked
    // is called the heartbeatInterval. The time step used should be either the heartbeat interval or the
    // expsosure time whatever is smallest. 

    double timeStep = min(telescope.getHeartbeatInterval(), exposureTime);


    // Later on we will have to convert from magnitudes to fluxes. Precompute a constant prefactor.
    // fluxOfV0Star is the photon flux [photons/s/m^2/nm] for a V=0 G2V-star.
    // Units of fluxFactor: [photons/s]
  

    const double hc = Constants::CLIGHT * Constants::HPLANCK * 1.e9;
    const double referenceFlux = fluxOfV0Star * hc / referenceWavelength; 
    const double fluxFactor = referenceFlux * binwidth * telescope.getLightCollectingArea();

    vector<double> wavePrefactors;
    wavePrefactors.reserve(binnumber);
    wavePrefactors.clear();

    vector<double> photonEnergies;
    photonEnergies.reserve(binnumber);
    photonEnergies.clear();

    for (int i=0; i<binnumber; i++)
    {
        double wavelength = lowerWavelength + i*binwidth + binwidth/2;
        double Ephoton = hc / wavelength;
        photonEnergies.push_back(Ephoton);
        double fluxFactorWave = transmissionEfficiencySpectral[i] * pow(referenceWavelength / wavelength, 5);
        wavePrefactors.push_back(fluxFactorWave);
    } 

    // Update the internal clock

    internalTime = startTime;

    // Take the flux of point sources (stars) into account.
    // Break up the exposure time in small intervals (hearbeat intervals) to track jitter while exposing.

    while (internalTime < startTime + exposureTime)
    {
        // Update the time-dependent parameters (if any) of some classes to to their value at the current time. 

        this->updateParameters(internalTime);
        telescope.updateParameters(internalTime);
        detector.updateParameters(internalTime);
        sky.updateParameters(internalTime);

        // Loop over all stars in the catalog, and add their flux to the subfield

        unsigned int NstarsInSubfield = 0;

        for (unsigned int n = 0; n < Nstars; n++)
        {
            // Compute the focal plane coordinates (in [mm]) of this particular star
           
            unsigned long starID;
            double RA, dec, Vmag;

            tie(starID, RA, dec, Vmag) = sky.getSelectedStar(n);
            double tempStar = sky.getSelectedStarTemp(n);
            
            double Xmm, Ymm;
            tie(Xmm, Ymm) = skyToFocalPlaneCoordinates(RA, dec);

            // If required, include field distortion

            if (includeFieldDistortion)
            {
                tie(Xmm, Ymm) = undistortedToDistortedFocalPlaneCoordinates(Xmm, Ymm);
            }

            // Compute the flux [photons] of this star
            // Photons are always an integer number, so round down.
            // To do this, calculate the spectrally correct transmissionEfficiency and flux fraction binwise

            double flux = 0.;

            for (int i=0; i<binnumber; i++)
            {
            double tempFactor = (exp(hc/(referenceWavelength*Constants::KBOLTZMANN*tempStar)) - 1) / (exp(photonEnergies[i]/(Constants::KBOLTZMANN*tempStar)) - 1);
            flux += floor(fluxFactor * wavePrefactors[i] * pow(10.0, -0.4 * Vmag) * timeStep * tempFactor) * QESpectral[i] / meanQE;
            } 

            // Let the detector add the flux to the appropriate pixel. 
            // Detector.flux() returns the pixel coordinates to which the flux was added.

            bool isInSubfield;
            double rowPix, colPix;    // subfield (not CCD) pixel coordinates

            tie(isInSubfield, rowPix, colPix) = detector.addFlux(Xmm, Ymm, flux);

            // If the star is indeed in the subfield, collect the following information to later write to HDF5
            //    1) average (Xmm, Ymm) coordinates of the star during the exposure                   [mm]
            //    2) average (row, col) pixel coordinates of the star on the CCD during the exposure  [pix]
            //    3) the total number of photons gathered of this star during the exposure            [photons]
            //    4) the total number of times that the star was in the subfield during the exposure
            //
            // Note: Due to jitter, the star can move in and out the subfield during the exposure

            if (isInSubfield)
            {
                NstarsInSubfield++;

                // If this is the first time we encounter this startTime, initialise the information

                if (detectedStarInfo.find(startTime) == detectedStarInfo.end())
                {
                    detectedStarInfo[startTime][starID] = {{Xmm, Ymm, rowPix, colPix, flux, 1.0}};
                }
                else
                {
                    // If this is the first time that we encounter this star ID associated with this startTime,
                    // initialise the information. If not, just update the info.

                    if (detectedStarInfo[startTime].find(starID) == detectedStarInfo[startTime].end())
                    {
                        detectedStarInfo[startTime][starID] = {{Xmm, Ymm, rowPix, colPix, flux, 1.0}};
                    }
                    else
                    {
                        detectedStarInfo[startTime][starID][0] += Xmm;      // Will be used to compute average Xmm during the exposure
                        detectedStarInfo[startTime][starID][1] += Ymm;      // Will be used to compute average Ymm during the exposure
                        detectedStarInfo[startTime][starID][2] += rowPix;   // Will be used to compute average pixel row during the exposure
                        detectedStarInfo[startTime][starID][3] += colPix;   // Will be used to compute average pixel column during the exposure
                        detectedStarInfo[startTime][starID][4] += flux;     // Total flux
                        detectedStarInfo[startTime][starID][5] += 1;        // # of times a star was on the subfield during an exposure 
                    }
                }
            }
        }

        Log.debug("Camera: at time " + to_string(internalTime) + ": incremented flux of " + to_string(NstarsInSubfield) + " stars in subfield");


        // Update the clock. Normally with 'timeStep', but if adding timeStep would overstep
        // the total exposure time, take the small rest time instead (but update the internal time first!).

        internalTime += timeStep;
        timeStep = min(timeStep, startTime + exposureTime - internalTime);
    }

    // Take the flux of the stellar background and the zodiacal light into account. Use one value for the entire subfield.
    // A negative value for the user given sky background value [phot/pix/s] signals that we should compute it ourselves.
    // Note: - The output of sky.zodiacalFlux() is in [J s^{-1} m^{-2} sr^{-1} m^{-1}]
    //       - As wavelength range we take the entire throughput band.
    //       - Photons are always an integer number, thus round down.
    //
    // The small sky background contribution during the readout with open shutter is not taken into account here
    // but in the function Detector::applyOpenShutterSmearing().

    totalSkyBackground = 0.0;


    if (userGivenSkyBackground < 0.0)
    {
        const double energyOfOnePhoton = Constants::CLIGHT * Constants::HPLANCK / (throughputLambdaC * 1.e-9);                // [J]
        const double lambda1 = (throughputLambdaC - throughputBandwidth/2.0) * 1.e-9;                                         // [m]
        const double lambda2 = (throughputLambdaC + throughputBandwidth/2.0) * 1.e-9;                                         // [m]
    
        const double zodiacalFlux = sky.zodiacalFlux(centerRA, centerDec, lambda1, lambda2)                                   // [phot/exposure]
                                    * (exposureTime + readoutTimeBeforeNextExposure) * transmissionEfficiency * telescope.getLightCollectingArea()
                                    * detector.getSolidAngleOfOnePixel(plateScale) / energyOfOnePhoton; 

        const double stellarBackgroundFlux = sky.stellarBackgroundFlux(centerRA, centerDec, lambda1, lambda2)                 // [phot/exposure]
                                             * (exposureTime + readoutTimeBeforeNextExposure) * transmissionEfficiency * telescope.getLightCollectingArea()
                                             * detector.getSolidAngleOfOnePixel(plateScale) / energyOfOnePhoton;      


        totalSkyBackground = floor(zodiacalFlux + stellarBackgroundFlux);
        detector.addFlux(totalSkyBackground);

        Log.debug("Camera: zodiacal flux level in subfield = " + to_string(zodiacalFlux) + " photons/pixel/exposure");
        Log.debug("Camera: stellar background flux level in subfield = " + to_string(stellarBackgroundFlux) + " photons/pixel/exposure");
    }
    else
    {
        totalSkyBackground = floor(userGivenSkyBackground * exposureTime * transmissionEfficiency);
        detector.addFlux(totalSkyBackground);

        Log.debug("Camera: user-given sky background flux over exposure= " + to_string(userGivenSkyBackground * exposureTime) + " photons/pixel/exposure");
    }

    // Save the sky background value that we added. [photons/pix/exposure]

    skyBackgroundValues.push_back(totalSkyBackground);

    // Save the transmissionEfficiency value that was calculated for this exposure.

    transmissionEfficiencyValues.push_back(transmissionEfficiency);

    return;
}

