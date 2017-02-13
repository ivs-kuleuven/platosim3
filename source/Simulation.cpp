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

Simulation::Simulation(string inputFilename, string outputFilename) {
	// Parse the configuration parameters file

	Log.info("Simulation: reading the input parameters file");

	ConfigurationParameters configParams(inputFilename);

	// Check if the output HDF5 filename already exists. If so, complain.

	if (fileExists(outputFilename)) {
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

	if (useJitterFromFile) {
		jitterGenerator = new JitterFromFile(configParams);
	} else {
		jitterGenerator = new JitterFromRedNoise(configParams);
	}

	// Depending on what the user requested, define the proper telescope thermo-elastic drift generator

	if (useDriftFromFile) {
		driftGenerator = new ThermoElasticDriftFromFile(configParams);
	} else {
		driftGenerator = new ThermoElasticDriftFromRedNoise(configParams);
	}

	// Initialise the spacecraft components

	platform = new Platform(configParams, hdf5File, *jitterGenerator);
	telescope = new Telescope(configParams, hdf5File, *platform,
			*driftGenerator);
	sky = new Sky(configParams);
	camera = new Camera(configParams, hdf5File, *telescope, *sky);
	detector = new Detector(configParams, hdf5File, *camera);

	// Write the input parameters to the output HDF5 file

	writeInputParametersToHDF5(configParams);

}

/**
 * \brief      Destructor, release memory of all spacecraft components
 */
Simulation::~Simulation() {
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

void Simulation::configure(ConfigurationParameters &configParams) {
	exposureTime = configParams.getDouble("ObservingParameters/ExposureTime");
	Nexposures = configParams.getInteger("ObservingParameters/NumExposures");
	useJitterFromFile = configParams.getBoolean("Platform/UseJitterFromFile");
	includeFieldDistortion = configParams.getBoolean(
			"Camera/IncludeFieldDistortion"); // do we want to do this or should this be asked to Camera?
	useDriftFromFile = configParams.getBoolean("Telescope/UseDriftFromFile");
}

/**
 * \brief      Loop over all exposures
 *
 * \param[in]  startTime  begin time of the very first exposure. Time is expressed in seconds in the rest of the code.
 */
void Simulation::run(double startTime) {
	// Update the internal clock

	currentTime = startTime;

	// Ensure that the proper PSF is set for the detector

	detector->setPsfForSubfield();

	// Loop over all exposures

	for (int n = 0; n < Nexposures; n++) {
		Log.info(
				"Simulation: Starting exposure " + to_string(n) + " at time "
						+ to_string(currentTime));

		currentTime = detector->takeExposure(currentTime, exposureTime);
	}

	writeStarCatalogToHDF5();
}

void Simulation::writeVersionInformationToHDF5() {
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
void Simulation::writeStarCatalogToHDF5() {
	Log.info(
			"Simulation: writing info on detected stars to HDF5 in /StarCatalog");

	set<unsigned int> allStarIDs = camera->getAllStarIDs();

	// For all detected stars, copy the equatorial sky coordinates and the magnitude
	// from the user-given star catalog to the output HDF5 file in a custom group.

	hdf5File.createGroup("/StarCatalog");

	const int Nstars = allStarIDs.size();
	vector<unsigned int> starIDs(Nstars); // set<> is not contiguous, vector<> is. Needed for HDF5.
	vector<double> RA(Nstars);
	vector<double> dec(Nstars);
	vector<double> Vmag(Nstars);
	vector<double> rowPix(Nstars);
	vector<double> colPix(Nstars);
	vector<double> xFPmm(Nstars);
	vector<double> yFPmm(Nstars);

	double raOpticalAxis, decOpticalAxis;
	tie(raOpticalAxis, decOpticalAxis) =
			telescope->getInitialPointingCoordinates();

	double focalPlaneAngle = telescope->getInitialFocalPlaneOrientation();

	double xFPrad, yFPrad;

	if (!allStarIDs.empty()) {
		int k = 0;
		for (auto starID : allStarIDs) {
			starIDs[k] = starID;
			tie(RA[k], dec[k]) = sky->getCoordinatesOfStarWithID(starID,
					Angle::degrees); // be careful, ra & dec returned in degrees!
			Vmag[k] = sky->getVmagnitudeOfStarWithID(starID);
			tie(xFPmm[k], yFPmm[k]) = camera->skyToFocalPlaneCoordinates(
					deg2rad(RA[k]), deg2rad(dec[k]), raOpticalAxis,
					decOpticalAxis, focalPlaneAngle);

			if (includeFieldDistortion) {
				tie(xFPmm[k], yFPmm[k]) =
						camera->undistortedToDistortedFocalPlaneCoordinates(
								xFPmm[k], yFPmm[k]);
			}

			tie(rowPix[k], colPix[k]) = detector->focalPlaneToPixelCoordinates(
					xFPmm[k], yFPmm[k]);
			k++;
		}

		hdf5File.writeArray("StarCatalog/", "starIDs", starIDs.data(),
				starIDs.size());
		hdf5File.writeArray("StarCatalog/", "RA", RA.data(), RA.size());
		hdf5File.writeArray("StarCatalog/", "Dec", dec.data(), dec.size());
		hdf5File.writeArray("StarCatalog/", "Vmag", Vmag.data(), Vmag.size());
		hdf5File.writeArray("StarCatalog/", "xFPmm", xFPmm.data(),
				xFPmm.size());
		hdf5File.writeArray("StarCatalog/", "yFPmm", yFPmm.data(),
				yFPmm.size());
		hdf5File.writeArray("StarCatalog/", "colPix", colPix.data(),
				colPix.size());
		hdf5File.writeArray("StarCatalog/", "rowPix", rowPix.data(),
				rowPix.size());
	} else {
		Log.warning(
				"Simulation: no information about detected stars to write to HDF5");
	}

}

/**
 * \brief [brief description]
 *
 */

void Simulation::writeInputParametersToHDF5(
		ConfigurationParameters &configParams) {
	Log.info("Simulation: writing input parameters to HDF5");

	// Make the parent group

	string parentGroup = "/InputParameters";
	hdf5File.createGroup(parentGroup);

	string subGroup;

	// Define some Lambda functions that will make it much easier to add the input parameters

	auto addDouble =
			[&] (string attributeName)
			{
				hdf5File.writeAttribute(parentGroup + "/" + subGroup, attributeName, configParams.getDouble(subGroup + "/" + attributeName));
			};

	auto addInteger =
			[&] (string attributeName)
			{
				hdf5File.writeAttribute(parentGroup + "/" + subGroup, attributeName, configParams.getInteger(subGroup + "/" + attributeName));
			};

	auto addLong =
			[&] (string attributeName)
			{
				hdf5File.writeAttribute(parentGroup + "/" + subGroup, attributeName, configParams.getLong(subGroup + "/" + attributeName));
			};

	auto addString =
			[&] (string attributeName)
			{
				hdf5File.writeAttribute(parentGroup + "/" + subGroup, attributeName, configParams.getString(subGroup + "/" + attributeName));
			};

	auto addBoolean =
			[&] (string attributeName)
			{
				hdf5File.writeAttribute(parentGroup + "/" + subGroup, attributeName, configParams.getBoolean(subGroup + "/" + attributeName));
			};

	auto addDoubleVector =
			[&] (string attributeName)
			{
				hdf5File.writeAttribute(parentGroup + "/" + subGroup, attributeName, configParams.getDoubleVector(subGroup + "/" + attributeName));
			};

	// Copy the input parameters to the output HDF5 file

	subGroup = "ObservingParameters";
	hdf5File.createGroup(parentGroup + "/" + subGroup);
	addInteger("NumExposures");
	addDouble("ExposureTime");
	addDouble("RApointing");
	addDouble("DecPointing");
	addDouble("Fluxm0");
	addDouble("SkyBackground");
	addString("StarCatalogFile");

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
	addDouble("AzimuthAngle");
	addDouble("TiltAngle");
	addDouble("LightCollectingArea");
	addDouble("TransmissionEfficiency");
	addBoolean("UseDrift");
	addBoolean("UseDriftFromFile");
	addDouble("DriftYawRms");
	addDouble("DriftPitchRms");
	addDouble("DriftRollRms");
	addDouble("DriftTimeScale");
	addString("DriftFileName");

	subGroup = "Camera";
	hdf5File.createGroup(parentGroup + "/" + subGroup);
	addDouble("FocalPlaneOrientation");
	addDouble("PlateScale");
	addDouble("FocalLength");
	addDouble("ThroughputBandwidth");
	addDouble("ThroughputLambdaC");
	addBoolean("IncludeFieldDistortion");
	subGroup = "Camera/FieldDistortion";
	hdf5File.createGroup(parentGroup + "/" + subGroup);
	addString("Type");
	addInteger("Degree");
	addDoubleVector("Coefficients");
	addDoubleVector("InverseCoefficients");

	subGroup = "PSF";
	hdf5File.createGroup(parentGroup + "/" + subGroup);
	addString("Model");
	subGroup = "PSF/Gaussian";
	hdf5File.createGroup(parentGroup + "/" + subGroup);
	addDouble("Sigma");
	addInteger("NumberOfPixels");

	subGroup = "PSF/FromFile";
	hdf5File.createGroup(parentGroup + "/" + subGroup);
	addString("Filename");
	addDouble("DistanceToOA");
	addDouble("RotationAngle");
	addInteger("NumberOfPixels");

	subGroup = "FEE";
	hdf5File.createGroup(parentGroup + "/" + subGroup);
	addDouble("NominalOperatingTemp");
	addDouble("ReadoutNoise");

	subGroup = "FEE/Gain";
	hdf5File.createGroup(parentGroup + "/" + subGroup);
	addDouble("RefValue");
	addDouble("Stability");
	addDouble("Delta");

	subGroup = "FEE/ElectronicOffset";
	hdf5File.createGroup(parentGroup + "/" + subGroup);
	addInteger("RefValue");
	addDouble("Stability");

	subGroup = "CCD";
	hdf5File.createGroup(parentGroup + "/" + subGroup);
	addDouble("OriginOffsetX");
	addDouble("OriginOffsetY");
	addDouble("Orientation");
	addInteger("NumColumns");
	addInteger("NumRows");
	addDouble("PixelSize");
	addLong("FullWellSaturation");
	addInteger("DigitalSaturation");
	addDouble("ReadoutNoise");
	addDouble("ReadoutTime");
	addDouble("FlatfieldPtPNoise");
	addDouble("NominalOperatingTemp");
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
	addBoolean("IncludeDigitalSaturation");
	addBoolean("WriteSubPixelImagesToHDF5");

	subGroup = "CCD/Gain";
	hdf5File.createGroup(parentGroup + "/" + subGroup);
	addDouble("RefValue");
	addDouble("Stability");
	addDouble("Delta");

	subGroup = "CCD/Vignetting";
	hdf5File.createGroup(parentGroup + "/" + subGroup);
	addDouble("ExpectedValue");

	subGroup = "CCD/QuantumEfficiency";
	hdf5File.createGroup(parentGroup + "/" + subGroup);
	addDouble("Efficiency");
	addDouble("RefAngle");
	addDouble("ExpectedValue");

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
	addLong("CTESeed");
	addLong("DriftSeed");
	addLong("FeeGainSeed");
	addLong("CcdGainSeed");
}
