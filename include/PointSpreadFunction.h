#ifndef PSF_H
#define PSF_H

#include <string>
#include <vector>

#include "Logger.h"
#include "HDF5File.h"
#include "ConfigurationParameters.h"

using namespace std;



class PointSpreadFunction
{
    public:

        PointSpreadFunction(ConfigurationParameters &);
        ~PointSpreadFunction();

    protected:
        void rotate(double angle);
        void configure(ConfigurationParameters &);


    private:
        void select();
        void rebin();

        // Determine if a psf has been selected
        bool isSelected = false;

        // The selected psf is copied into this array
        arma::Mat<float> psfMap;

        // The dimensions of the selected psf map will be set here
        unsigned int numRowsPsfMap;
        unsigned int numColumnsPsfMap;

        // The HDF5 file that holds the PSFs
        HDF5File *hdf5file;

        // Loaded from the configuration, i.e. PSF/Filename
        string location;

        // Name of the HDF5 group that contains the PSF datasets
        string groupName;
};



#endif