#ifndef PSF_H
#define PSF_H

#include <string>
#include <vector>
#include <algorithm>

#include "armadillo"

#include "ArrayOperations.h"
#include "ConfigurationParameters.h"
#include "Exceptions.h"
#include "FileUtilities.h"
#include "StringUtilities.h"
#include "HDF5File.h"
#include "Logger.h"
#include "Units.h"

using namespace std;





class PointSpreadFunction
{
    public:

        PointSpreadFunction(ConfigurationParameters &);
        ~PointSpreadFunction();

        void rotate(double angle);
        void select(double radius);
        void rebin(unsigned int numSubPixelsPerPixel);

        arma::Mat<float> getPsfMap();

    protected:
        void configure(ConfigurationParameters &);


    private:

        // Determine if a psf has been selected
        bool isSelected = false;

        // Determine if this psf has been rotated
        bool isRotated = false;

        // Determine if this psf has been rotated
        bool isRebinned = false;

        // The angle by which the PSF is rotated with respect to the positive x-axis.
        // Positive angles rotated counter-clockwise.
        double rotationAngle = 0;    // [radians]

        // The selected psf is copied into this array
        arma::Mat<float> psfMap;

        // The HDF5 file that holds the PSFs
        HDF5File hdf5file;

        // Loaded from the configuration, i.e. PSF/Filename
        string absolutePath;

        // Name of the HDF5 group that contains the PSF datasets
        string groupName;

        // Number of sub-pixels per pixel that was used to generate the PSF
        unsigned int numberOfSubPixelsPerPixel;

        // Location of the reference subpixel, i.e. the center of the PSF
        double xCenter;
        double yCenter;
};






namespace psfdata
{
    const arma::vec radius = {
        0.0000, 1.4141, 2.8273, 4.2388, 5.6477, 7.0532, 8.4545, 9.8508, 11.2413, 12.6253, 14.0019, 15.3707, 16.7308, 18.0817, 18.8876
    };
}

#endif