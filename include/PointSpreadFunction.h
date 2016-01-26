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


    private:
        void loadConfiguration(ConfigurationParameters &);
        void select();
        void rotate();
        void rebin();

        // The HDF5 file that holds the PSFs
        HDF5File *hdf5file;

        // Loaded from the configuration, i.e. PSFFileName
        string location;

        // Name of the HDF5 group that contains the PSF datasets
        string groupName;
};



#endif