/**
 * \class PointSpreadFunction
 * 
 * \brief Base class for the PSF, both symmetrical and asymmetrical. 
 */

#include "PointSpreadFunction.h"





/**
 * \brief Default constructor.
 * 
 * \param configParam: Configuration parameters.
 * 
 * \param hdf5file: HDF5 file from which to read the PSF.
 */
PointSpreadFunction::PointSpreadFunction(ConfigurationParameters &configParam, HDF5File &hdf5file) : HDF5Writer(hdf5file){}





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







/*
 * This function gets the selected psf
 */

arma::fmat PointSpreadFunction::getOriginalPSF()
{
  psfMap = ArrayOperations::rotateArray(psfMap, -rotationAngle);
  psfMap /= arma::accu(psfMap);
  return psfMap;
}



