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
: HDF5Writer(hdf5file) {}








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
