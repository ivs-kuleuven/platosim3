/**
 * \class     PointSpreadFunction 
 * 
 * \brief     The PointSpreadFunction provides the PSF at different positions in the field.
 * 
 * \details
 * 
 * The PSF intensity maps are provided as a HDF5 file which contains several groups of PSFs for 
 * different temperature and different positions in the field. Currently only PSFs for 6000K are
 * available. 
 * 
 * 
 * \todo
 * 
 * There are two hardcoded values used in this class, i.e. the top-level groupName of the HDF5 PSF file, 
 * and the name of the dataset containing the PSF array. The latter should actually be generated from the 
 * ra, dec from the center of the sub-field, but this information is currently not passed into the select()
 * method.
 * 
 */

#include "PointSpreadFunction.h"








/**
 * \brief      Constructor
 *
 * \details
 * 
 * Load the configuration parameters and initialize the internal variables. The default group 
 * for which PSF will be loaded is currently set to "6000" (which is the only group available at the moment).
 * 
 * The HDF5 PSF file is loaded and some basic checking is done.
 * 
 * \param      configParam  configuration parameters for the PSF
 */

PointSpreadFunction::PointSpreadFunction(ConfigurationParameters &configParam, HDF5File &hdf5file)
: HDF5Writer(hdf5file)
{
    // Parse the parameters from the configuration file.

    configure(configParam);

    // Create the groups in the HDF5 file where the different PSFs and their description will be saved. - Moved after configure for spectral dependency
    initHDF5Groups();

    isSelected = false;
    isRotated = false;
    
    if (isLoadedFromFile)
    {
        // Prepare the psfFile by performing some basic checks

        if ( ! FileUtilities::fileExists(absolutePath) )
        {
            throw FileException("PointSpreadFunction: trying to load the PSF HDF5 file (" + absolutePath + "), but file doesn't exist.");
        }
    
        try
        {
            psfFile.open(absolutePath);
        }
        catch(H5::FileIException ex)
        {
            Log.error("H5::FileIException: " + string(ex.getCDetailMsg()));
            throw H5FileException("PointSpreadFunction: Could not open HDF5 file: " + absolutePath);
        }
    
        string groupName = "wave-1";  //%% Spectral dependency, wave-1 is the combined PSF of all wavelengths

        if ( ! psfFile.hasGroup(groupName) )
        {
            throw H5FileException("PointSpreadFunction: The HDF5 file (" + absolutePath + ") doesn't contain the expected group \"" + groupName + "\".");
        }
    
        Log.info("PointSpreadFunction: Opened the HDF5 file " + absolutePath + " containing the PSFs");
    }

}








/**
 * \brief      Destructor
 * 
 * \details
 * 
 * Close the HDF5 file and release the memory.
 * 
 */
PointSpreadFunction::~PointSpreadFunction()
{
    flushOutput();
    psfFile.close();
}









/**
 * \brief Write all recorded and remaining information to the HDF5 output file.
 */
void PointSpreadFunction::flushOutput()
{
    Log.info("PointSpreadFunction: Flushing output to HDf5 file.");

    if (! isSelected)
        return;
    
    // Save the PSF subpixel map when it is rebinned to pixel level.

   // rebinToPixels();  //%% Needs addtional work for the subsubfields
}









/**
 * \brief Creates the group(s) in the HDF5 file where the PSF information will be stored. 
 *        These group(s) have to be created once, at the very beginning.
 */
void PointSpreadFunction::initHDF5Groups()
{
    Log.debug("PointSpreadFunction: initialising HDF5 groups");

    hdf5File.createGroup("/PSF");

    for (int binnumber=0; binnumber<wave_bins; binnumber++)  //%% Spectral dependency: create a group for each wavelength
    {
	string group;
	if (wave_bins == 1)  //%% If only one wavelength is chosen use the multicolor PSF
	{
	group = "/PSF/wavebin" + to_string(binnumber-1);
	}
	else
	{
	group = "/PSF/wavebin" + to_string(binnumber);
	}
	hdf5File.createGroup(group);
    }
}









/**
 * \brief      Return a two-dimensional Gaussian function
 *
 * \details    The Gaussian PSF is calculated from the following equation:
 *
 *             \f[
 *                f(x, y)  =  \frac{1}{2\pi\sigma^{2}} e^{-[(x - \mu_x)^{2} + (y - \mu_y)^{2}] / (2\sigma^{2})}
 *             \f]
 *            
 *             where \f$\mu_x\f$ and \f$\mu_y\f$ are the mean (center points) and 
 *             \f$\sigma\f$ is the standard deviation.
 *            
 *             The Gaussian PSF is centered in an array at subpixel level, i.e. the array
 *             has dimension [numberOfSubPixelsPerPixel * numberOfPixels,
 *             numberOfSubPixelsPerPixel * numberOfPixels], but the standard deviation 
 *             \f$\sigma\f$ is defined at pixel level.
 *
 * \return     a 2D Gaussian PSF
 */
arma::fmat PointSpreadFunction::getGaussianPsf()
{

    double centerRow = numberOfSubPixelsPerPixel * numberOfPixels / 2.0;
    double centerColumn = numberOfSubPixelsPerPixel * numberOfPixels / 2.0;

    if (sigma < 1.0 / numberOfSubPixelsPerPixel)
    {
        throw IllegalArgumentException("The width of the Gaussian PSF must be larger than the size of a subpixel.");
    }

    // Generate the Gaussian PSF at the sub-pixel level

    double width = sigma * numberOfSubPixelsPerPixel;
    double denominator = 2.0 * width * width;
    double sumPSF = 0.0;

    arma::fmat gaussianPsf (numberOfPixels * numberOfSubPixelsPerPixel, numberOfPixels * numberOfSubPixelsPerPixel, arma::fill::zeros);

    for (unsigned int xi = 0; xi < gaussianPsf.n_cols; xi++)
    {
        for (unsigned int xj = 0; xj < gaussianPsf.n_rows; xj++)
        {
            gaussianPsf(xj, xi) = exp( - (pow(xi - centerColumn, 2) + pow(xj - centerRow, 2)) / denominator);
            sumPSF += gaussianPsf(xj, xi);
        }
    }

    // Normalize the gaussian, so that the flux is conserved
    
    gaussianPsf /= sumPSF;

    return gaussianPsf;
}










/**
 * \brief      Configure the PointSpreadFunction object using the
 *             ConfigurationParameters
 *
 * \param[in]  configParam  the configuration parameters
 */
void PointSpreadFunction::configure(ConfigurationParameters &configParam)
{
    string model = configParam.getString("PSF/Model");

    wave_bins = configParam.getInteger("Camera/WavelengthBins");  //%% Spectral dependence: Read how many wavelength bins are to be processed

    if (model == "MappedGaussian")
    {
        // The user specified to use a Gaussian shape PSF
        // The number of sub-pixels per pixel that will be used to calculate the 
        // Gaussian PSF is equal to the number of sub-pixels per pixels for the sub-field.

        isGaussian                = true;
        sigma                     = configParam.getDouble("PSF/MappedGaussian/Sigma");
        numberOfPixels            = configParam.getInteger("PSF/MappedGaussian/NumberOfPixels");

        // The Gaussian PSF shall be created with a resolution equal to that of the sub-field
        
        numberOfSubPixelsPerPixel = configParam.getInteger("SubField/SubPixels");
    } 
    else if (model == "MappedFromFile")
    {
        // The user specified to use the pre-calculated PSFs from file
        // The number of sub-pixels per pixel is derived from the number of pixels specified 
        // and the size of the array in the file. (The number of sub-pixels should be in the file).
      
        isLoadedFromFile          = true;
        absolutePath              = configParam.getAbsoluteFilename("PSF/MappedFromFile/Filename");
        numberOfPixels            = configParam.getInteger("PSF/MappedFromFile/NumberOfPixels");
        requestedDistanceToOA     = deg2rad(configParam.getDouble("PSF/MappedFromFile/DistanceToOA"));
        requestedRotationAngle    = deg2rad(configParam.getDouble("PSF/MappedFromFile/RotationAngle"));

    numsubsubfieldsx 	   = configParam.getInteger("SubField/NumSubSubFieldsRows");  //%% Multiple fields, for spectral dependence
    numsubsubfieldsy 	   = configParam.getInteger("SubField/NumSubSubFieldsColumns");  //% Multiple fields, for spectral dependence

    }
    else
    {
        string errorMessage = "PointSpreadFunction: Model '" + model + "' is not supported.";
        Log.error(errorMessage);
        throw IllegalArgumentException(errorMessage);
    }


}








/**
 * \brief      Return the distance to the optical axis as requested by the user.
 * 
 * \details    When the returned value is negative [-1], the user input will be ignored
 *
 * \return     distance to optical axis [rad]
 */
double PointSpreadFunction::getRequestedDistanceToOpticalAxis()
{
    return requestedDistanceToOA;
}









/**
 * \brief      Return the orientation angle of the PSF as requested by the user.
 *
 * \details    This orientation is what the user specified in the input file and
 *             is not necessarily equal to the rotation angle of the PSF (which
 *             is also compensated for the CCD orientation).
 *
 * \return     the orientation as requested by the user
 */
double PointSpreadFunction::getRequestedRotationAngle()
{
    return requestedRotationAngle;
}









/**
 * \brief      Select the proper PSF matching the given radius closest.
 *
 * \details    The generated PSF is position dependent and a look-up table has been provided that 
 *             contains the field radial coordinate along the line of sight for which the PSF was generated.
 *  
 *             The appropriate PSF will be selected, i.e. the PSF for which the angular distance of 
 *             the centre of the sub-field from the centre of the focal plane matches best.
 * 
 * \param[in]  radius  angular separation of the source for which to select the PSF [radians]
 */
void PointSpreadFunction::select(double radius, int fieldnumber, int fieldmax)
{
    using StringUtilities::dtos;
    
    // TODO: Should we take any action if different PSFs are selected for this object?

    if (isSelected)
    {
        Log.warning("PointSpreadFunction: Another PSF was previously selected.");

psfVector.clear();  //%% Clear the vectors as we use append
rotationVector.clear();
    }

    if (isGaussian)
    {
        psfMap = getGaussianPsf();
        isSelected = true;
        rotationAngle = 0.0;

        hdf5File.writeAttribute("/PSF", "selectedPSF", "Gaussian PSF selected with sigma=" + to_string(sigma));

        return;
    }

    radius = rad2deg(radius);

//%% Read radii from psf hdf5, added for spectral dependency
    arma::vec radii;
    int field = 1;  //%% Starting radii at 1 so use 0 as test case
    while (field != 0)
    {
	string groupName = "wave-" +to_string(1) + "/radius" + to_string(field);  //%% Take the radii from the multicolor PSF which has to be there
	if (psfFile.hasGroup(groupName))
	{
	    double rad = psfFile.readDoubleGroupAttribute(groupName, "radius[deg]");
	    int sz = radii.size();
	    radii.resize(sz+1);
	    radii(sz) = rad;
	    field = field+1;
	}
	else
	{
	    field = 0;
	}
    }

    // Convert radius into the string angularRadiusGroup that identifies the psf dataset in the HDF5 file.
    // We work with a lookup table psfdata::radius which contains fixed radius values for which PSF data
    // was generated. The algorithm is to select the PSF with radius closest to the given radius by 
    // subtracting the given radius from the tabulated radius and then selecting the lowest value to find the index.
    
    arma::vec rads = radii - radius;  //%% Use radii, containing all radii read from the hdf5

    rads = abs(rads);

    arma::uword index;
    rads.min(index);

    if (index > radii.n_elem-1)  //%% radii vector for spectral dependence
    {
        Log.warning("PointSpreadFunction: Radius index (" + to_string(index) + ") is out of bounds.");
        index = radii.n_elem-1;  //%% radii vector for spectral dependence
    }

    stringstream myStream;

    const string angularRadiusGroup = "radius" + to_string(index+1);  //%% Changed string for spectral dependency

    const string azimuthDataset = "0";  // TODO: hardcoded value!

//%% Spectral dependency: Loop over all wavebins
    for (int binnumber=0; binnumber<wave_bins; binnumber++)
{

    string temperatureGroup;
    if (wave_bins == 1){temperatureGroup = "wave" + to_string(binnumber-1);}  //%% Use multicolor if only one wavebin
    else{temperatureGroup = "wave" + to_string(binnumber);}  //% Otherwise use correct temperature, curretntly only works if identical to binpositions

    string groupName = temperatureGroup + "/" + angularRadiusGroup;

    if ( ! psfFile.hasGroup(groupName) )
    {
        throw FileException("PointSpreadFunction: The HDF5 file (" + absolutePath + ") doesn't contain the expected group \"" + groupName + "\".");
    }

    if ( ! psfFile.hasDataset(groupName, azimuthDataset) )
    {
        throw FileException("PointSpreadFunction: The HDF5 file (" + absolutePath + ") doesn't contain the expected dataset \"" + azimuthDataset + "\".");
    }

    // Load the psf array into the psfMap
    
    psfFile.readArray("/" + groupName, azimuthDataset, psfMap);

    psfVector.push_back(psfMap);  //%% Save the psf in a vector to keep all wavelengths

    // The PSFs that are currently used are rotated with respect to the focal plane x-axis.
    // The rotation angle is given as an attribute to the dataset that contains the PSF.

    double angle = psfFile.readDoubleDatasetAttribute(groupName, azimuthDataset, "orientation");

    rotationAngle = deg2rad(angle);

    rotationVector.push_back(rotationAngle);  //%% keep the rotation angles of all bins
    string group;
    if (wave_bins == 1)
    {
	group = "/PSF/wavebin" + to_string(binnumber-1);
    }
    else
    {
	group = "/PSF/wavebin" + to_string(binnumber);
    }
    Log.debug("PointSpreadFunction: Selected PSF " + groupName + "/" + azimuthDataset + " for wavebin " + to_string(binnumber) + ", rotation set to " + dtos(angle) + " degrees.");
    hdf5File.writeAttribute(group, "selectedPSFfield" + to_string(fieldnumber), "Realistic PSF selected from group " + groupName + "/" + azimuthDataset + ".");
}

    // We should be able to read the number of sub-pixels per pixel that was used to generate the PSFs
    // from an attribute in the HDF5 file. Unfortunately, this is not available and we therefore derive 
    // the number from the array size of the psfMap and the number of pixels, currently specified 
    // in the input file.

    numberOfSubPixelsPerPixel = psfMap.n_rows / numberOfPixels;

    if (fieldnumber == fieldmax){  //%% If processed the last subsubfield, then all psfs have been chosen, for spectral dependency
      isSelected = true;
    }
}







/**
 * \brief      Rotate the PSF with the given angle.
 * 
 * \details
 * 
 * Beware that the PSF that has been provided is already rotated, this will be taken into account.
 *
 * \param[in]  angle  angle by which the PSF should be rotated [radians]
 */
void PointSpreadFunction::rotate(double angle, int fieldnumber, int fieldmax)
{

    // We do not need to rotate a Gaussian PSF
    // Even if we do not rotate the Gaussian PSF, we do save the psfMap as a rotatedPSF.
    // This is to keep consistency in the output file where we do not save the selected PSF,
    // but we do save the rotated PSF.

    if (isGaussian)
    {
        hdf5File.writeArray("/PSF", "rotatedPSF", psfMap);
        hdf5File.writeAttribute("/PSF", "rotationAngle", rotationAngle);
        return;
    }

    if (isRotated)
    {
        Log.warning("PointSpreadFunction: Ignoring rotation: PSF was already rotated before.");
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
	Log.debug("PointSpreadFunction: rotated current PSF over angle " + to_string(rad2deg(newAngle)) + " deg for wavebin " + to_string(binnumber));
}

    if (fieldnumber == fieldmax){
	isRotated = true;  //%% If processed the last subsubfield, then all psfs have been rotated, for spectral dependency
    }

        // Write the psfMap of the rotated PSF to the HDF5 output file  //%% Not done yet

//        hdf5File.writeArray("/PSF", "rotatedPSF", psfMap);
//        hdf5File.writeAttribute("/PSF", "rotationAngle", rotationAngle);

    }
}











/**
 * @brief      Rebin the PSF map to the target number of subpixels
 *
 * @details    The number of subpixels used to generate the PSF is not
 *             necessarily the same as the number of subpixels per pixel for the
 *             detector. So, the PSF needs to be rebinned to the number of
 *             subpixels per pixel for the detector, which is given as
 *             targetSubPixels.
 *
 * @note       This method does not change the psfMap of the PointSpreadFunction class.
 * 
 * @param[in]  targetSubPixels  the target number of subpixels
 *
 * @return     the rebinned PSF map
 */
vector<arma::Mat<float>> PointSpreadFunction::rebinToSubPixels(unsigned int targetSubPixels, int fieldnumber)  //%% Changed to vector for spectral dependency
{
    if (targetSubPixels == numberOfSubPixelsPerPixel)
	return psfVector;    //%% Changed to vector for spectral dependency

//%% Added a loop over all wavelength bins
    vector<arma::Mat<float>> rebinnedVector;
    arma::fmat rebinnedMap;
    for (int binnumber=0; binnumber<wave_bins; binnumber++)
{
    unsigned int binSize = psfVector[binnumber].n_rows / numberOfSubPixelsPerPixel * targetSubPixels;
    rebinnedMap = ArrayOperations::rebin(psfVector[binnumber], binSize, binSize);  //%% Overwrite rebinned map

    rebinnedMap /= arma::accu(rebinnedMap);

    rebinnedVector.push_back(rebinnedMap);  //%% Keep the value of the bin in a vector

    string group;
    if (wave_bins == 1)
    {
	group = "/PSF/wavebin" + to_string(binnumber-1);
    }
    else
    {
	group = "/PSF/wavebin" + to_string(binnumber);
    }
    hdf5File.writeArray(group, "rebinnedPSFsubPixelfield" + to_string(fieldnumber), rebinnedMap);
}
    isRebinned = true;

    // Write the rebinned PSF to the output HDF5 file  //%% Not done yet

//    hdf5File.writeArray("/PSF", "rebinnedPSFsubPixel", rebinnedMap);

    return rebinnedVector;  //%% Changed to vector for spectral dependency
}









/**
 * @brief      Rebin the PSF map to the target number of pixels
 *
 * @details    The PSF subpixel map will be rebinned to a pixel map.
 * 
 * @note       This method does not change the psfMap of the PointSpreadFunction class.
 *
 * @param[in]  targetPixels  the target number of pixels
 * 
 * @return     the rebinned PSF map
 * 
 */

vector<arma::fmat> PointSpreadFunction::rebinToPixels()  //%% Changed to vector for spectral dependency
{
    vector<arma::fmat> rebinnedVector;
    for (int binnumber=0; binnumber<wave_bins; binnumber++)
{
    string group;
    if (wave_bins == 1)
    {
	group = "/PSF/wavebin" + to_string(binnumber-1);
    }
    else
    {
	group = "/PSF/wavebin" + to_string(binnumber);
    }

    unsigned int binSize = psfVector[binnumber].n_rows / numberOfSubPixelsPerPixel;

    arma::fmat rebinnedMap = ArrayOperations::rebin(psfVector[binnumber], binSize, binSize);

    rebinnedMap /= arma::accu(rebinnedMap);

    rebinnedVector.push_back(rebinnedMap);

    // Write the rebinned PSF to the output HDF5 file //%% Not Done yet

    //hdf5File.writeArray(group, "rebinnedPSFpixel", rebinnedMap);
}
    isRebinned = true;

    return rebinnedVector;
}
