/**
 * \class PointSpreadFunction
 *
 * \brief Base class for the PSF, both symmetrical and asymmetrical.
 */

#include "PointSpreadFunction.h"





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
PointSpreadFunction::PointSpreadFunction(ConfigurationParameters &configParam, HDF5File &hdf5file) : HDF5Writer(hdf5file)
{
    // Create the groups in the HDF5 file where the different PSFs and their description will be saved.

    initHDF5Groups();

    // Parse the parameters from the configuration file.

    configure(configParam);

    isSelected = false;
    isRotated = false;

    // Prepare the psfFile by performing some basic checks

    if (!FileUtilities::fileExists(absolutePath))
    {
        throw FileException("PointSpreadFunction: trying to load the PSF HDF5 file (" + absolutePath + "), but file doesn't exist.");
    }

    try
    {
        psfFile.open(absolutePath);
    }
    catch (H5::FileIException ex)
    {
        Log.error("H5::FileIException: " + string(ex.getCDetailMsg()));
        throw H5FileException("PointSpreadFunction: Could not open HDF5 file: " + absolutePath);
    }

    Log.info("PointSpreadFunction: Opened the HDF5 file " + absolutePath + " containing the PSFs");

    initializeDistortionMap();
}





/**
 * \brief Destructor.
 *
 * \details Closes the HDF5 file and releases the memory.
 */
PointSpreadFunction::~PointSpreadFunction()
{
    flushOutput();
    psfFile.close();
}





/**
 * \brief Configures the PointSpreadFunction object using the given
 *        configuration parameters
 *
 * \param[in] configParam: Configuration parameters.
 */
void PointSpreadFunction::configure(ConfigurationParameters &configParam)
{
    string model = configParam.getString("PSF/Model");

    if (model == "MappedFromFile")
    {
        // The user specified to use the pre-calculated PSFs from file
        // The number of sub-pixels per pixel is derived from the number of pixels specified
        // and the size of the array in the file. (The number of sub-pixels should be in the file).

        absolutePath = configParam.getAbsoluteFilename("PSF/MappedFromFile/Filename");
        numberOfPixels = configParam.getInteger("PSF/MappedFromFile/NumberOfPixels");
	writeHighResolutionPSF = configParam.getBoolean("ControlHDF5Content/WriteHighResolutionPSF");
    }
    else
    {
        string errorMessage = "PointSpreadFunction: Model '" + model + "' is not supported.";
        Log.error(errorMessage);
        throw IllegalArgumentException(errorMessage);
	writeHighResolutionPSF = false;
    }

    string focalLengthSource = configParam.getString("Camera/FocalLength/Source");
    if (focalLengthSource == "ConstantValue")
    {
        focalLengthValue = configParam.getDouble("Camera/FocalLength/ConstantValue") * 1000;     // [m] -> [mm]
    }
    else if (focalLengthSource == "FromFile")
    {
        string focalLengthInputFile = configParam.getAbsoluteFilename("Camera/FocalLength/FromFile");
        Parameter<double> *focalLength = new Parameter<double>(focalLengthInputFile, 1000);             // [m] -> [mm]
        focalLengthValue = (*focalLength)();
    }
    Log.info("PointSpreadFunction: Using a focal length: " + to_string(focalLengthValue) + " mm");
}



/**
 * \brief Creates the group(s) in the HDF5 file where the PSF information will be stored.
 *        These group(s) has (have) to be created once, at the very beginning.
 */
void PointSpreadFunction::initHDF5Groups()
{
    Log.debug("PointSpreadFunction: initialising HDF5 groups");

    hdf5File.createGroup("/PSF");
}





/**
 * \brief Writes all recorded and remaining information to the HDF5 output file.
 */
void PointSpreadFunction::flushOutput()
{
    Log.info("PointSpreadFunction: Flushing output to HDf5 file.");

    if (!isSelected)
        return;

    // Save the PSF subpixel map when it is rebinned to pixel level.

    rebinToPixels();
}




/**
 * \brief Selects the proper PSF matching the given focal-plane coordinates closest.  The appropriate PSF will be
 *        selected, i.e. the PSF for which the focal-plane position matches best.
 *
 * \param[in] xFP: Focal-plane x-coordinate [mm].
 *
 * \param[in] yFP: Focal-plane y-coordinate [mm].
 */
void PointSpreadFunction::select(double xFP, double yFP)
{
    using StringUtilities::dtos;

    // TODO: Should we take any action if different PSFs are selected for this object?

    if (isSelected)
    {
        Log.warning("PointSpreadFunction: Another PSF was previously selected.");
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
            // The following two attributes are stored as strings in the HDF5 file rather than doubles,
            // so we need to convert them.

            xPsf = stod(psfFile.readStringDatasetAttribute("/", datasetName, "centerCoordinates1"));
            yPsf = stod(psfFile.readStringDatasetAttribute("/", datasetName, "centerCoordinates2"));

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



    rotationAngle = 0.0;

    // We should be able to read the number of sub-pixels per pixel that was used to generate the PSFs
    // from an attribute in the HDF5 file. Unfortunately, this is not available and we therefore derive
    // the number from the array size of the psfMap and the number of pixels, currently specified
    // in the input file.

    numberOfSubPixelsPerPixel = psfMap.n_rows / numberOfPixels;

    hdf5File.writeAttribute("/PSF", "selectedPSF", "Realistic PSF selected from dataset " + datasetName + ".");

    isSelected = true;
}






/**
 * \brief Rotates the PSF with the given angle.  Beware that the PSF that has been provided is already
 *        rotated, this will be taken into account.
 *
 * \param[in] angle: Angle by which the PSF should be rotated [radians].
 */
void PointSpreadFunction::rotate(double angle)
{
    if (isRotated)
    {
        Log.warning("PointSpreadFunction: Ignoring rotation: PSF was already rotated before.");
    }
    else
    {
        // This part of the code currently assumes that the PSF has not been rotated yet.
        // The rotationAngle of the generated PSFs is 45 degrees, so the actual rotation should be
        // the requested angle minus the rotationAngle of the PSF.

        psfMap = ArrayOperations::rotateArray(psfMap, angle);

        rotationAngle = angle;
        isRotated = true;

        psfMap /= arma::accu(psfMap);

        Log.debug("PointSpreadFunction: rotated current PSF over angle " + to_string(rad2deg(angle)) + " deg");

        // Write the psfMap of the rotated PSF to the HDF5 output file
	if (writeHighResolutionPSF)
	{
        hdf5File.writeArray("/PSF", "highResPSF", psfMap);
        hdf5File.writeAttribute("/PSF", "rotationAngle", rotationAngle);
	}
    }
}


/** \brief Rebins the PSF sub-pixel map to pixel level.
 *        This method does not change the psfMap of the PointSpreadFunction class.
 *
 * \param[in] targetPixels: Target number of pixels.
 *
 * \return Rebinned PSF map.
 *
 */
arma::fmat PointSpreadFunction::rebinToPixels()
{
    unsigned int binSize = psfMap.n_rows / numberOfSubPixelsPerPixel;

    arma::fmat rebinnedMap = ArrayOperations::rebin(psfMap, binSize, binSize);

    isRebinned = true;

    rebinnedMap /= arma::accu(rebinnedMap);

    // Write the rebinned PSF to the output HDF5 file

    hdf5File.writeArray("/PSF", "rebinnedPSFpixel", rebinnedMap);

    return rebinnedMap;
}





/**
 * \brief Rebins the PSF map to the target number of sub-pixels. The number of sub-pixels used
 *        to generate the PSF is not necessarily the same as the number of sub-pixels per pixel
 *        for the detector. So the PSF needs to be rebinned to the number of sub-pixels per pixel
 *        for the detector, which is the specified target number of sub-pixels.
 *        This method does not change the psfMap of the PointSpreadFunction class.
 *
 * \param[in] targetSubPixels: Target number of sub-pixels (i.e. after rebinning).
 *
 * \return Rebinned PSF map (with the given number of sub-pixels).
 */
arma::fmat PointSpreadFunction::rebinToSubPixels(unsigned int targetSubPixels)
{
    if (targetSubPixels == numberOfSubPixelsPerPixel)
        return psfMap;

    unsigned int binSize = psfMap.n_rows / numberOfSubPixelsPerPixel * targetSubPixels;

    arma::fmat rebinnedMap = ArrayOperations::rebin(psfMap, binSize, binSize);

    isRebinned = true;

    rebinnedMap /= arma::accu(rebinnedMap);

    // Write the rebinned PSF to the output HDF5 file

    hdf5File.writeArray("/PSF", "rebinnedPSFsubPixel", rebinnedMap);

    return rebinnedMap;
}







/**
 * This function gets the selected psf
 */

arma::fmat PointSpreadFunction::getOriginalPSF()
{
  psfMap = ArrayOperations::rotateArray(psfMap, -rotationAngle);
  psfMap /= arma::accu(psfMap);
  return psfMap;
}












/**
 * Initialize the distortion map. This can be done either by copying it from
 * the HDF5 file or generated if the starPointing angles are given.
 */

void PointSpreadFunction::initializeDistortionMap()
{

    // First we check if a Coordinates Map is defined in the hdf5 file

    if (psfFile.hasGroup("/Coordinates map"))
    {
        readDistortionmapFromFile();
    }
    else
    {
        // If no map is defined we will generate such a map.
        generateDistortionMap();
    }
}








/**
 * \brief Return a vector containing the coordinates for the best fitting
 *        Wang distortion model.
 *
 * \param[in] focalLength: The focalLength that we use.
 *
 * \return vector<double> with the 6 coefficient to descibe the Wang model.
 */
vector<double> PointSpreadFunction::estimateDistortionCoefficients(double focalLength)
{
    MappedDistortion distortion = MappedDistortion(xFP, yFP, xFPdist, yFPdist, focalLength);
    vector<double> coefficients = distortion.getParameters();

    return coefficients;
}







/**
 * \brief Return a vector containing the coordinates for the best fitting
 *        Wang inverse distortion model.
 *
 * \param[in] focalLength: The focalLength that we use.
 *
 * \return vector<double> with the 6 coefficient to descibe the Wang model for
 *         the inverse distortion.
 */
vector<double> PointSpreadFunction::estimateInverseDistortionCoefficients(double focalLength)
{
    MappedDistortion distortion = MappedDistortion(xFPdist, yFPdist, xFP, yFP, focalLength);
    std::vector<double> coefficients = distortion.getParameters();

    return coefficients;
}











/**
 * \brief Return a vector containing the coordinates for the best fitting
 *        5th order 2D polynomial distortion.
 *
 * \return vector<vector<double>> containing two vectors, each containing the
 *         36 coordinates to describe the distortion model.
 */
vector<vector<double>> PointSpreadFunction::estimatePolynomialCoefficients()
{

    MappedDistortion distortion = MappedDistortion(xFP, yFP, xFPdist, yFPdist);

    vector<double> coefficientsX = distortion.getParametersX();
    vector<double> coefficientsY = distortion.getParametersY();

    vector<vector<double>> coefficients = {coefficientsX, coefficientsY};
    return coefficients;
}





/**
 * \brief Return a vector containing the coordinates for the best fitting
 *        5th order 2D polynomial inverse distortion.
 *
 * \return vector<vector<double>> containing two vectors, each containing the
 *         36 coordinates to describe the inverse distortion model.
 */
vector<vector<double>> PointSpreadFunction::estimateInversePolynomialCoefficients()
{

    MappedDistortion distortion = MappedDistortion(xFPdist, yFPdist, xFP, yFP);

    vector<double> coefficientsX = distortion.getParametersX();
    vector<double> coefficientsY = distortion.getParametersY();

    vector<vector<double>> coefficients = {coefficientsX, coefficientsY};
    return coefficients;
}










/**
 * Reads the distortion map from the HDf5 and copies its content in distortionMap.
 */
void PointSpreadFunction::readDistortionmapFromFile()
{
    // We read in the table that converts the undistorted coordinates to distorted coordinates from the psf HDF5 file.

    psfFile.readArray("/Coordinates map/Undistorted", "x", xFP);
    psfFile.readArray("/Coordinates map/Undistorted", "y", yFP);
    psfFile.readArray("/Coordinates map/Distorted", "x", xFPdist);
    psfFile.readArray("/Coordinates map/Distorted", "y", yFPdist);

}






/**
 * Generates a distortion map if no map is given in the HDf5 file.
 */
void PointSpreadFunction::generateDistortionMap()
{

    unsigned int datasetIndex = 1;
    string datasetName;
    string selectedDatasetName;

    while(true)
    {
        datasetName = to_string(datasetIndex);

        if(psfFile.hasDataset("/", datasetName))
        {
            // The distorted FP coordinats are stored as strings in the HDF5 file rather than doubles,
            // so we need to convert them.

            double xPsf = stod(psfFile.readStringDatasetAttribute("/", datasetName, "centerCoordinates1"));
            double yPsf = stod(psfFile.readStringDatasetAttribute("/", datasetName, "centerCoordinates2"));

            // We then read in the star pointing angles and from those we obtain the undistorted FP coordiantes.

            double starPointing1[1];
            double starPointing2[1];
            psfFile.readArrayDatasetAttribute("/", datasetName, "starPointing1", starPointing1);
            psfFile.readArrayDatasetAttribute("/", datasetName, "starPointing2", starPointing2);

            double xPsfUndistorted = focalLengthValue*tan(deg2rad(starPointing1[0]));
            double yPsfUndistorted = focalLengthValue*tan(deg2rad(starPointing2[0]));

            // Once we have these values we add them to their corresponding vectors
            const std::array<double, 4> coordinates = { xPsfUndistorted, yPsfUndistorted, xPsf, yPsf };
            xFP.push_back(xPsfUndistorted);
            yFP.push_back(yPsfUndistorted);
            xFPdist.push_back(xPsf);
            yFPdist.push_back(yPsf);
        }
        else
        {
            break;
        }

        datasetIndex++;
    }
}
