
#ifndef HDF5WRITER_H
#define HDF5WRITER_H

#include <limits>
#include <string>
#include "hdf5file.h"

using namespace std;



class HDF5Writer
{
    
    public:

        HDF5Writer(HDF5File &hdf5File, string groupName);
        ~HDF5Writer();

        virtual void flushOutput();

    protected:


    private:

        HDF5File &hdf5File;
        string hdf5GroupName;
    
};

#endif
