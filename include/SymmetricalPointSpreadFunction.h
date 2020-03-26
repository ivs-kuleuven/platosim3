#ifndef SYMMETRICALPOINTSPREADFUNCTION_H
#define SYMMETRICALPOINTSPREADFUNCTION_H

#include <string>
#include <vector>
#include <algorithm>

#include "armadillo"

#include "ConfigurationParameters.h"
#include "Constants.h"
#include "Exceptions.h"
#include "FileUtilities.h"
#include "StringUtilities.h"
#include "Logger.h"
#include "Units.h"
#include "PointSpreadFunction.h"

using namespace std;





class SymmetricalPointSpreadFunction : public PointSpreadFunction
{
    public:

        SymmetricalPointSpreadFunction(ConfigurationParameters &configParam, HDF5File &hdf5File);

        void select(double radius, int fieldnumber, int fieldmax);

        double getRequestedDistanceToOpticalAxis();
        double getRequestedRotationAngle();

        arma::fmat getGaussianPsf();

        virtual void rotate(double angle, int fieldnumber, int fieldmax);

    protected:

        virtual void configure(ConfigurationParameters &configParam);
 
    private:

       
        // Indicates whether the PSF is Gaussian, because this is treated differently
        // (not location dependent, no rotation needed)
        
        bool isGaussian = false;

        // The PSF shall be loaded from an HDF5 file
        // This option can not be true if isGaussian is already true!

        bool isLoadedFromFile = false;
        
        // Name of the HDF5 group that contains the PSF datasets

        string groupName;

        // Width (standard deviation) of the Gaussian PSF [pixels]

        double sigma;

        // Angular distance to the optical axis as requested by the user.  A negative value 
        // indicates no user input, i.e. auto-compute.

        double requestedDistanceToOA = -1.0;     // [radians]

        // Angle by which the PSF should be rotated with respect to the positive x-axis.
        // Positive angles rotate counter-clockwise.  A negative value indicates no user input,
        // i.e. auto-compute

        double requestedRotationAngle = -1.0;    // [radians]
};



#endif
