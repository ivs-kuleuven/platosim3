/**
 * \class AsymmetricalPointSpreadFunction 
 * 
 * \brief Class for rotationally asymmetrical PSFs, at different location on the focal plane.
 * 
 * \details The PSF intensity maps are provided as an HDF5 file, which contains one group per location on the
 *          focal plane.
 */
#include "AsymmetricalPointSpreadFunction.h"





/**
 * \brief Constructor.
 * 
 * \details Initialises the groups in the HDF5 file where the different maps (i.e. pixel map,
 *          bias register map, smearing map, etc.) will be saved. 
 * 
 * The following maps are initialized to zero (partly through the base class Detector):
 *      - pixelMap 
 *      - subPixelMap
 *      - biasMap
 *      - smearingMap
 *      - flatfieldMap
 *      - throughputMap
 *      - cteMap
 * 
 * The flatfieldMap is filled at sub-pixel level, the throughputMap and cteMap are filled at pixel level.
 * 
 * \param configParam: Configuration parameters for the detector.
 * 
 * \param hdf5file:HFD5 file to write the detector images to.
 * 
 * \param camera: Camera to which to attach the detector.
 * 
 * \param readoutTimeBeforeNextExposure Duration of the readout that takes place before the next exposure can start [s].
 */
AsymmetricalPointSpreadFunction::AsymmetricalPointSpreadFunction(ConfigurationParameters &configParam, HDF5File &hdf5file) : PointSpreadFunction(configParam, hdf5file)
{
    // Parse the parameters from the configuration file.

    configure(configParam);

    // Create the groups in the HDF5 file where the different PSFs and their description will be saved.

    initHDF5Groups();

    isSelected = false;
    isRotated = false;

    // Prepare the psfFile by performing some basic checks

    if (!FileUtilities::fileExists(absolutePath))
    {
        throw FileException("AsymmetricalPointSpreadFunction: trying to load the PSF HDF5 file (" + absolutePath + "), but file doesn't exist.");
    }

    try
    {
        psfFile.open(absolutePath);
    }
    catch (H5::FileIException ex)
    {
        Log.error("H5::FileIException: " + string(ex.getCDetailMsg()));
        throw H5FileException("AsymmetricalPointSpreadFunction: Could not open HDF5 file: " + absolutePath);
    }

    Log.info("AsymmetricalPointSpreadFunction: Opened the HDF5 file " + absolutePath + " containing the PSFs");
}





/**
 * \brief Configures the AsymmetricalPointSpreadFunction object using the given
 *        configuration parameters
 *
 * \param[in] configParam: Configuration parameters.
 */
void AsymmetricalPointSpreadFunction::configure(ConfigurationParameters &configParam)
{
    string model = configParam.getString("PSF/Model");

    wave_bins = configParam.getInteger("Camera/WavelengthBins");  //%% Spectral dependence: Read how many wavelength bins are to be processed

    if (model == "MappedFromFileAsymmetrical")
    {
        // The user specified to use the pre-calculated PSFs from file
        // The number of sub-pixels per pixel is derived from the number of pixels specified 
        // and the size of the array in the file. (The number of sub-pixels should be in the file).

        absolutePath = configParam.getAbsoluteFilename("PSF/MappedFromFileAsymmetrical/Filename");
        numberOfPixels = configParam.getInteger("PSF/MappedFromFileAsymmetrical/NumberOfPixels");

    numsubsubfieldsx 	   = configParam.getInteger("SubField/NumSubSubFieldsRows");  //%% Multiple fields, for spectral dependence
    numsubsubfieldsy 	   = configParam.getInteger("SubField/NumSubSubFieldsColumns");  //% Multiple fields, for spectral dependence
    }
    else
    {
        string errorMessage = "AsymmetricalPointSpreadFunction: Model '" + model + "' is not supported.";
        Log.error(errorMessage);
        throw IllegalArgumentException(errorMessage);
    }
}





/**
 * \brief Selects the proper PSF matching the given focal-plane coordinates closest.  The appropriate PSF will be 
 *        selected, i.e. the PSF for which the focal-plane position matches best.
 * 
 * \param[in] xFP: Focal-plane x-coordinate [mm].
 * 
 * \param[in] yFP: Focal-plane y-coordinate [mm].
 */
void AsymmetricalPointSpreadFunction::select(double xFP, double yFP, int fieldnumber, int fieldmax)
{
    using StringUtilities::dtos;
    
    // TODO: Should we take any action if different PSFs are selected for this object?

    if (isSelected)
    {
        Log.warning("PointSpreadFunction: Another PSF was previously selected.");
psfVector.clear();  //%% Clear the vectors as we use append
rotationVector.clear();
    }

    unsigned int datasetIndex = 1;
    string datasetName;
    double xPsf, yPsf, distanceSquared;
    double minDistanceSquared = std::numeric_limits<double>::max();
    string selectedDatasetName;

    while(true)
    {
        datasetName = to_string(datasetIndex);

        if(psfFile.hasDataset("/", datasetName))
        {
            xPsf = psfFile.readDoubleDatasetAttribute("/", datasetName, "centerCoordinates1");
            yPsf = psfFile.readDoubleDatasetAttribute("/", datasetName, "centerCoordinates2");

            distanceSquared = pow(xPsf - xFP, 2) + pow(yPsf - yFP, 2);

            if(distanceSquared < minDistanceSquared)
            {
                minDistanceSquared = distanceSquared;
                selectedDatasetName = datasetName;
            }
        }
        else
        {
            break;
        }

        datasetIndex++;
    }

    psfFile.readArray("/", selectedDatasetName, psfMap);

for (int binnumber=0; binnumber<wave_bins; binnumber++)
{
psfVector.push_back(psfMap);  //%% Save the psf in a vector to keep all wavelengths, quick and dirty fix

    rotationAngle = 0.0;

rotationVector.push_back(rotationAngle);  //%% keep the rotation angles of all bins, quick and dirty fix
}

    // We should be able to read the number of sub-pixels per pixel that was used to generate the PSFs
    // from an attribute in the HDF5 file. Unfortunately, this is not available and we therefore derive 
    // the number from the array size of the psfMap and the number of pixels, currently specified 
    // in the input file.

    numberOfSubPixelsPerPixel = psfMap.n_rows / numberOfPixels;

    hdf5File.writeAttribute("/PSF", "selectedPSF", "Realistic PSF selected from dataset " + datasetName + ".");

    if (fieldnumber == fieldmax){  //%% If processed the last subsubfield, then all psfs have been chosen, for spectral dependency
      isSelected = true;
    }
}





/** 
 * \brief Rotates the PSF with the given angle.  Beware that the PSF that has been provided is already 
 *        rotated, this will be taken into account.
 *
 * \param[in] angle: Angle by which the PSF should be rotated [radians].
 */
void AsymmetricalPointSpreadFunction::rotate(double angle, int fieldnumber, int fieldmax)
{
    if (isRotated)
    {
        Log.warning("AsymmetricalPointSpreadFunction: Ignoring rotation: PSF was already rotated before.");
    }
    else
    {
        // This part of the code currently assumes that the PSF has not been rotated yet.
        // The rotationAngle of the generated PSFs is 45 degrees, so the actual rotation should be
        // the requested angle minus the rotationAngle of the PSF.

//%% Added loop for spectral dependence, rotate all bins
        for (int binnumber=0; binnumber<wave_bins; binnumber++)
{
	double newAngle = angle - rotationVector[binnumber];

	psfVector[binnumber] = ArrayOperations::rotateArray(psfVector[binnumber], newAngle);  //%% Changed to vector for spectral dependency

        rotationAngle = newAngle;  

	psfVector[binnumber] /= arma::accu(psfVector[binnumber]);		//NORMALIZING  //%% Changed to vector for spectral dependency

    string group;
    if (wave_bins == 1)
    {
	group = "/PSF/wavebin" + to_string(binnumber-1);
    }
    else
    {
	group = "/PSF/wavebin" + to_string(binnumber);
    } 
        hdf5File.writeArray(group, "rotatedPSFfield" + to_string(fieldnumber), psfVector[binnumber]);
        hdf5File.writeAttribute(group, "rotationAnglefield" + to_string(fieldnumber), rotationAngle);
	Log.debug("AsymmetricalPointSpreadFunction: rotated current PSF over angle " + to_string(rad2deg(newAngle)) + " deg for wavebin " + to_string(binnumber));
}

    if (fieldnumber == fieldmax){
	isRotated = true;  //%% If processed the last subsubfield, then all psfs have been rotated, for spectral dependency
    }

        // Write the psfMap of the rotated PSF to the HDF5 output file  //%% Not done yet

//        hdf5File.writeArray("/PSF", "rotatedPSF", psfMap);
//        hdf5File.writeAttribute("/PSF", "rotationAngle", rotationAngle);

    }
}
