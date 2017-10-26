/**
 * \class Simulation
 * 
 * \brief The starting point for any simulation.
 * 
 */

#include "Simulation.h"


/**
 * \brief      Constructor
 * 
 * \details
 * 
 * The constructor reads the YAML input file, and creates the HDF5 output file.
 * Based on the user input a Jitter generator is created and all spacecraft
 * components are initialized.
 *
 * \param[in]  inputFilename   the YAML input file 
 * \param[in]  outputFilename  the HDF5 output file
 */

Simulation::Simulation(string inputFilename, string outputFilename)
{
    // Parse the configuration parameters file

    Log.info("Simulation: reading the input parameters file");

    ConfigurationParameters configParams(inputFilename);

    // Check if the output HDF5 filename already exists. If so, complain.

    if (fileExists(outputFilename))
    {
        Log.error("Simulation: Output file name already exists. Aborting.");
        exit(1);
    }

    // Open the HDF5 output file where the images will be written

    hdf5File.open(outputFilename);

    // Write the version info to the output HDF5 file

    writeVersionInformationToHDF5();

    // Configure the Simulation object using the configuration parameters file

    configure(configParams);

    // Depending on what the user requested, define the proper platform jitter generator

    if (useJitterFromFile)
    {
        jitterGenerator = new JitterFromFile(configParams);
    }
    else
    {
        jitterGenerator = new JitterFromRedNoise(configParams);
    }

    // Depending on what the user requested, define the proper telescope thermo-elastic drift generator

    if (useDriftFromFile)
    {
        driftGenerator = new ThermoElasticDriftFromFile(configParams);
    }
    else
    {
        driftGenerator = new ThermoElasticDriftFromRedNoise(configParams);
    }

    if(useFeeTemperatureFromFile)
    {
    		feeTemperatureGenerator = new TemperatureFromFile(configParams, "FEE");
    }
    else if(useFeeNominalTemperature)
    {
    		feeTemperatureGenerator = new NominalTemperature(configParams, "FEE");
    }

    if(useDetectorTemperatureFromFile)
    {
    		detectorTemperatureGenerator = new TemperatureFromFile(configParams, "CCD");
    }
    else if(useDetectorNominalTemperature)
    {
    		detectorTemperatureGenerator = new NominalTemperature(configParams, "CCD");
    }

    // Initialise the spacecraft components

    platform   = new Platform(configParams, hdf5File, *jitterGenerator);
    telescope  = new Telescope(configParams, hdf5File, *platform, *driftGenerator);
    sky        = new Sky(configParams);
    camera     = new Camera(configParams, hdf5File, *platform, *telescope, *sky);

    // Depending on how the PSF is computed (analytically or pre-mapped) the Detector object is different.

    if ((psfModel == "MappedGaussian") || (psfModel == "MappedFromFile"))
    {
        detector = new DetectorWithMappedPSF(configParams, hdf5File, *camera, *feeTemperatureGenerator, *detectorTemperatureGenerator);
    }
    else if (psfModel == "AnalyticGaussian")
    {
        detector = new DetectorWithAnalyticGaussianPSF(configParams, hdf5File, *camera, *feeTemperatureGenerator, *detectorTemperatureGenerator);
    }
    else if (psfModel == "AnalyticNonGaussian")
    {
        detector = new DetectorWithAnalyticNonGaussianPSF(configParams, hdf5File, *camera, *feeTemperatureGenerator, *detectorTemperatureGenerator);
    }
    else
    {
        string errorMessage = "Simulation: PSF Model '" + psfModel + "' is not supported.";
        Log.error(errorMessage);
        throw IllegalArgumentException(errorMessage);
    }

    // Write the input parameters to the output HDF5 file

    writeInputParametersToHDF5(configParams);

}







/**
 * \brief      Destructor, release memory of all spacecraft components
 */
Simulation::~Simulation()
{
    // Delete order is the inverse of the order in which they were created

    delete detector;
    delete camera;
    delete telescope;
    delete sky;
    delete platform;
    delete jitterGenerator;
    delete driftGenerator;
    
    // Close the output hdf5 file

    hdf5File.close();
}








/**
 * \brief Configure the Simulation object using the input parameter file
 * 
 * \param configParams  Contains all configuration parameters from the input file
 */

void Simulation::configure(ConfigurationParameters &configParams)
{
    exposureTime      = configParams.getDouble("ObservingParameters/ExposureTime"); 
    beginExposureNr   = configParams.getInteger("ObservingParameters/BeginExposureNr");
    numExposures      = configParams.getInteger("ObservingParameters/NumExposures");
    useJitterFromFile = configParams.getBoolean("Platform/UseJitterFromFile");
    includeFieldDistortion = configParams.getBoolean("Camera/IncludeFieldDistortion"); // do we want to do this or should this be asked to Camera?
    useDriftFromFile  = configParams.getBoolean("Telescope/UseDriftFromFile");  
    psfModel          = configParams.getString("PSF/Model");
    useFeeTemperatureFromFile = configParams.getString("FEE/Temperature") == "FromFile";
    useFeeNominalTemperature = configParams.getString("FEE/Temperature") == "Nominal";
    useDetectorTemperatureFromFile = configParams.getString("CCD/Temperature") == "FromFile";
    useDetectorNominalTemperature = configParams.getString("CCD/Temperature") == "Nominal";
    readoutTime       = configParams.getDouble("CCD/ReadoutTime"); 
}








/**
 * \brief      Loop over all exposures
 *
 * \param[in]  startTime  begin time of the very first exposure. Time is expressed in seconds in the rest of the code.
 */

void Simulation::run()
{
    // Update the internal clock

    currentTime = beginExposureNr * (exposureTime + readoutTime);

    // Loop over all exposures

    for (int n = beginExposureNr; n < beginExposureNr + numExposures; n++)
    {
        Log.info("Simulation: Starting exposure " + to_string(n) + " at time " + to_string(currentTime) );
        
        currentTime = detector->takeExposure(n, currentTime, exposureTime);
    }

    writeStarCatalogToHDF5();
}










/**
 * \brief Take care that the version of the simulator is included in the HDF5 file,.
 */

void Simulation::writeVersionInformationToHDF5()
{
    Log.info("Simulation: writing version information to HDF5");

    // Make the parent group

    string parentGroup = "/Version";
    hdf5File.createGroup(parentGroup);
 
    hdf5File.writeAttribute(parentGroup, "Application", string("PlatoSim3"));
    hdf5File.writeAttribute(parentGroup, "GitVersion", string(GIT_DESCRIBE));

}









/**
 * \brief      Write information about the stars that were detected in the subField 
 *             to the HDF5 output file. 
 * 
 * \details    The Camera collects all the stars that fall within the boundaries of the subField.
 * 
 *             This function should only be called after all exposures have been taken in order 
 *             to have the complete collections of stars that have been detected in the subField.
 *             
 */
void Simulation::writeStarCatalogToHDF5()
{
    Log.info("Simulation: writing info on detected stars to HDF5 in /StarCatalog");

    set<unsigned int> allStarIDs = camera->getAllStarIDs();

    // For all detected stars, copy the equatorial sky coordinates and the magnitude 
    // from the user-given star catalog to the output HDF5 file in a custom group.
    
    hdf5File.createGroup("/StarCatalog");

    const int Nstars = allStarIDs.size();
    vector<unsigned int> starIDs(Nstars);    // set<> is not contiguous, vector<> is. Needed for HDF5.
    vector<double> RA(Nstars);
    vector<double> dec(Nstars);
    vector<double> Vmag(Nstars);
    vector<double> rowPix(Nstars);
    vector<double> colPix(Nstars);
    vector<double> xFPmm(Nstars);
    vector<double> yFPmm(Nstars);

    double xFPrad, yFPrad;

    if (!allStarIDs.empty())
    {
        int k = 0;
        for (auto starID: allStarIDs)
        {
            starIDs[k] = starID;
            tie(RA[k], dec[k]) = sky->getCoordinatesOfStarWithID(starID, Angle::degrees);  // be careful, ra & dec returned in degrees!
            Vmag[k] = sky->getVmagnitudeOfStarWithID(starID);
            const bool useInitialOrientation = true;
            tie(xFPmm[k], yFPmm[k]) = camera->skyToFocalPlaneCoordinates(deg2rad(RA[k]), deg2rad(dec[k]), useInitialOrientation);
            
            if (includeFieldDistortion)
            {
               tie(xFPmm[k], yFPmm[k]) = camera->undistortedToDistortedFocalPlaneCoordinates(xFPmm[k], yFPmm[k]);
            }

            tie(rowPix[k], colPix[k]) = detector->focalPlaneToPixelCoordinates(xFPmm[k], yFPmm[k]);
            k++;
        }

        hdf5File.writeArray("StarCatalog/", "starIDs", starIDs.data(), starIDs.size());
        hdf5File.writeArray("StarCatalog/", "RA",      RA.data(), RA.size());
        hdf5File.writeArray("StarCatalog/", "Dec",     dec.data(), dec.size());
        hdf5File.writeArray("StarCatalog/", "Vmag",    Vmag.data(), Vmag.size());
        hdf5File.writeArray("StarCatalog/", "xFPmm",    xFPmm.data(), xFPmm.size());
        hdf5File.writeArray("StarCatalog/", "yFPmm",    yFPmm.data(), yFPmm.size());
        hdf5File.writeArray("StarCatalog/", "colPix",    colPix.data(), colPix.size());
        hdf5File.writeArray("StarCatalog/", "rowPix",    rowPix.data(), rowPix.size());
    }
    else
    {
        Log.warning("Simulation: no information about detected stars to write to HDF5");
    }
}











/**
 * \brief [brief description]
 *
 */

void Simulation::writeInputParametersToHDF5(ConfigurationParameters &configParams)
{
    Log.info("Simulation: writing input parameters to HDF5");

    // Make the parent group

    string parentGroup = "/InputParameters";
    hdf5File.createGroup(parentGroup);

    string subGroup;

    // Define some Lambda functions that will make it much easier to add the input parameters

    auto addDouble = [&] (string attributeName) 
    { 
        hdf5File.writeAttribute(parentGroup + "/" + subGroup, attributeName, configParams.getDouble(subGroup + "/" + attributeName));
    };

    auto addInteger = [&] (string attributeName) 
    { 
        hdf5File.writeAttribute(parentGroup + "/" + subGroup, attributeName, configParams.getInteger(subGroup + "/" + attributeName));
    };

    auto addLong = [&] (string attributeName) 
    { 
        hdf5File.writeAttribute(parentGroup + "/" + subGroup, attributeName, configParams.getLong(subGroup + "/" + attributeName));
    };

    auto addString = [&] (string attributeName) 
    { 
        hdf5File.writeAttribute(parentGroup + "/" + subGroup, attributeName, configParams.getString(subGroup + "/" + attributeName));
    };

    auto addBoolean = [&] (string attributeName) 
    { 
        hdf5File.writeAttribute(parentGroup + "/" + subGroup, attributeName, configParams.getBoolean(subGroup + "/" + attributeName));
    };

    auto addDoubleVector = [&] (string attributeName) 
    {
        hdf5File.writeAttribute(parentGroup + "/" + subGroup, attributeName, configParams.getDoubleVector(subGroup + "/" + attributeName));
    };

    auto addIntegerVector = [&] (string attributeName) 
    {
        hdf5File.writeAttribute(parentGroup + "/" + subGroup, attributeName, configParams.getIntegerVector(subGroup + "/" + attributeName));
    };


    // Copy the input parameters to the output HDF5 file

    subGroup = "ObservingParameters";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addDouble("MissionDuration");
    addInteger("NumExposures");
    addInteger("BeginExposureNr");
    addDouble("ExposureTime");
    addDouble("RApointing");
    addDouble("DecPointing");
    addDouble("Fluxm0");
    addString("StarCatalogFile");

    subGroup = "Sky";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addDouble("SkyBackground");
    addBoolean("IncludeCosmics");
    subGroup = "Sky/Cosmics";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addDouble("CosmicHitRate");
    addDoubleVector("CosmicTrailLength");
    addDoubleVector("CosmicIntensity");

    subGroup = "Platform";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addBoolean("UseJitter");
    addBoolean("UseJitterFromFile");
    addDouble("JitterYawRms");
    addDouble("JitterPitchRms");
    addDouble("JitterRollRms");
    addDouble("JitterTimeScale");
    addString("JitterFileName");

    subGroup = "Telescope";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addString("GroupID");
    addDouble("AzimuthAngle");
    addDouble("TiltAngle");
    addDouble("LightCollectingArea");
    addBoolean("UseDrift");
    addBoolean("UseDriftFromFile");
    addDouble("DriftYawRms");
    addDouble("DriftPitchRms");
    addDouble("DriftRollRms");
    addDouble("DriftTimeScale");
    addString("DriftFileName");
    subGroup = "Telescope/TransmissionEfficiency";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addDouble("BOL");
    addDouble("EOL");

    subGroup = "Camera";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addDouble("PlateScale");
    addDouble("ThroughputBandwidth");
    addDouble("ThroughputLambdaC");
    addBoolean("IncludeFieldDistortion");
    subGroup = "Camera/FieldDistortion";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addString("Type");
    addString("Source");
    addDoubleVector("ConstantCoefficients");
    addDoubleVector("ConstantInverseCoefficients");
    addString("CoefficientsFromFile");
    addString("InverseCoefficientsFromFile");
    subGroup = "Camera/FocalPlaneOrientation";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addString("Source");
    addDouble("ConstantValue");
    addString("FromFile");
    subGroup = "Camera/FocalLength";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addString("Source");
    addDouble("ConstantValue");
    addString("FromFile");

    subGroup = "PSF";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addString("Model");
    subGroup = "PSF/MappedGaussian";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addDouble("Sigma");
    addInteger("NumberOfPixels");

    subGroup = "PSF/MappedFromFile";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addString("Filename");
    addDouble("DistanceToOA");
    addDouble("RotationAngle");
    addInteger("NumberOfPixels");

    subGroup = "PSF/AnalyticGaussian";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addDouble("Sigma00");
    addDouble("SigmaX18");
    addDouble("SigmaY18");

    subGroup = "PSF/AnalyticNonGaussian";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addString("ParameterFileName");
    subGroup = "PSF/AnalyticNonGaussian/Sigma";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addString("Source");
    addDouble("ConstantValue");
    addString("FromFile");

	subGroup = "FEE";
	hdf5File.createGroup(parentGroup + "/" + subGroup);
	addDouble("NominalOperatingTemperature");
	addString("Temperature");
	addString("TemperatureFileName");
	addDouble("ReadoutNoise");

	subGroup = "FEE/Gain";
	hdf5File.createGroup(parentGroup + "/" + subGroup);
	addDouble("RefValue");
	addDouble("Stability");
	addDouble("ThreeSigma");

	subGroup = "FEE/ElectronicOffset";
	hdf5File.createGroup(parentGroup + "/" + subGroup);
	addInteger("RefValue");
	addDouble("Stability");

	subGroup = "CCD";
	hdf5File.createGroup(parentGroup + "/" + subGroup);
    addString("Position");
	addDouble("OriginOffsetX");
	addDouble("OriginOffsetY");
	addDouble("Orientation");
	addInteger("NumColumns");
    addInteger("NumRows");
    addInteger("FirstRowExposed");
	addDouble("PixelSize");
	addLong("FullWellSaturation");
	addInteger("DigitalSaturation");
	addDouble("ReadoutNoise");
    addDouble("ReadoutTime");
    addDouble("FlatfieldPtPNoise");
	addDouble("NominalOperatingTemperature");
	addString("Temperature");
	addString("TemperatureFileName");
    addBoolean("IncludeFlatfield");
    addBoolean("IncludePhotonNoise");
    addBoolean("IncludeReadoutNoise");
    addBoolean("IncludeCTIeffects"); 
    addBoolean("IncludeOpenShutterSmearing");
    addBoolean("IncludeVignetting");
    addBoolean("IncludePolarization");
    addBoolean("IncludeParticulateContamination");
    addBoolean("IncludeMolecularContamination");
    addBoolean("IncludeQuantumEfficiency");
    addBoolean("IncludeConvolution");
    addBoolean("IncludeFullWellSaturation");
    addBoolean("IncludeQuantisation");
    addBoolean("IncludeDigitalSaturation");
    addBoolean("WriteSubPixelImagesToHDF5");

	subGroup = "CCD/Gain";
	hdf5File.createGroup(parentGroup + "/" + subGroup);
	addDouble("RefValue");
	addDouble("Stability");
	addDouble("ThreeSigma");

    subGroup = "CCD/Vignetting";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addDouble("ExpectedValue");

    subGroup = "CCD/QuantumEfficiency";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
//    addDouble("Efficiency");
    addDouble("RefAngle");
    addDouble("RelativeRefEfficiency");
    addDouble("MeanQuantumEfficiency");
//    addDouble("ExpectedValue");

    subGroup = "CCD/Polarization";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addDouble("Efficiency");
    addDouble("RefAngle");
    addDouble("ExpectedValue");

    subGroup = "CCD/Contamination";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addDouble("ParticulateContaminationEfficiency");
    addDouble("MolecularContaminationEfficiency");

    subGroup = "CCD/CTI";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addString("Model");
    subGroup = "CCD/CTI/Simple";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addDouble("MeanCTE");
    subGroup = "CCD/CTI/Short2013";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addDouble("Beta");
    addDouble("Temperature");
    addInteger("NumTrapSpecies");
    addDoubleVector("TrapDensity");
    addDoubleVector("TrapCaptureCrossSection");
    addDoubleVector("ReleaseTime");

    subGroup = "SubField";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addInteger("ZeroPointRow");
    addInteger("ZeroPointColumn");
    addInteger("NumColumns");
    addInteger("NumRows");
    addInteger("NumBiasPrescanRows");
    addInteger("NumSmearingOverscanRows");
    addInteger("SubPixels");

    subGroup = "RandomSeeds";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addLong("ReadOutNoiseSeed");
    addLong("PhotonNoiseSeed");
    addLong("JitterSeed");
    addLong("FlatFieldSeed");
    addLong("DriftSeed");
	addLong("FeeGainSeed");
	addLong("CcdGainSeed");

    subGroup = "CameraGroups";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addDoubleVector("AzimuthAngle");
    addDoubleVector("TiltAngle");

    subGroup = "CCDPositions";
    hdf5File.createGroup(parentGroup + "/" + subGroup);
    addDoubleVector("OriginOffsetX");
    addDoubleVector("OriginOffsetY");
    addDoubleVector("Orientation");
    addIntegerVector("NumColumns");
    addIntegerVector("NumRows");
    addIntegerVector("FirstRowForNormalCamera");
    addIntegerVector("FirstRowForFastCamera");


}
