#ifndef PSF_H
#define PSF_H

#include <string>
#include <vector>
#include <algorithm>

#include "armadillo"

#include "ArrayOperations.h"
#include "ConfigurationParameters.h"
#include "Constants.h"
#include "Exceptions.h"
#include "FileUtilities.h"
#include "StringUtilities.h"
#include "HDF5File.h"
#include "HDF5Writer.h"
#include "Logger.h"
#include "Units.h"

using namespace std;





class PointSpreadFunction : public HDF5Writer
{
    public:

        PointSpreadFunction(ConfigurationParameters &, HDF5File &hdf5File);
        virtual ~PointSpreadFunction();

        virtual void initHDF5Groups() override;
        virtual void flushOutput() override;

        void rotate(double angle);
        void select(double radius);

        double getRequestedDistanceToOpticalAxis();
        double getRequestedRotationAngle();

        arma::fmat rebinToSubPixels(unsigned int targetSubPixels);

    protected:
        void configure(ConfigurationParameters &);
        arma::fmat rebinToPixels();
        arma::fmat getGaussianPsf();

    private:

        // Determine if a psf has been selected
        bool isSelected = false;

        // Determine if this psf has been rotated
        bool isRotated = false;

        // Determine if this psf has been rotated
        bool isRebinned = false;

        // Handle a Gaussian PSF slightly different, i.e. not location dependent, no rotation needed
        bool isGaussian = false;

        // The PSF shall be loaded from an HDF5 file
        // This option can not be true if isGaussian is already true!
        bool isLoadedFromFile = false;
        
        // The selected psf is copied into this array
        arma::Mat<float> psfMap;

        // The HDF5 file that holds the PSFs
        HDF5File psfFile;

        // Loaded from the configuration, i.e. PSF/Filename
        string absolutePath;

        // Name of the HDF5 group that contains the PSF datasets
        string groupName;

        // Number of sub-pixels per pixel that was/is used to generate the PSF
        unsigned int numberOfSubPixelsPerPixel;

        // Number of pixels in the field that holds the PSF
        unsigned int numberOfPixels;

        // Width (standard deviation) of the Gaussian PSF [pixels]
        double sigma;

        // The actual rotation angle of the PSF with respect to the x-axis orientation of the focal plane
        double rotationAngle = 0.0;    // [radians]

        // The angular distance to the Optical Axis as requested by the user
        // A negative value indicates no user input,, i.e. auto-compute
        double requestedDistanceToOA = -1.0;     // [radians]

        // The angle by which the PSF should be rotated with respect to the positive x-axis.
        // Positive angles rotate counter-clockwise.
        // A negative value indicates no user input,, i.e. auto-compute
        double requestedRotationAngle = -1.0;    // [radians]

};






namespace psfdata
{
    const arma::vec radius = {
        0.0000, 1.4141, 2.8273, 4.2388, 5.6477, 7.0532, 8.4545, 9.8508, 11.2413, 12.6253, 14.0019, 15.3707, 16.7308, 18.0817, 18.8876
    };
}

#endif