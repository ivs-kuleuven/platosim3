
#ifndef HDF5WRITER_H
#define HDF5WRITER_H

#include <limits>
#include <string>
#include "hdf5file.h"

using namespace std;



class HDF5Writer
{
    
    public:

        HDF5Writer(HDF5File &hdf5File);
        ~HDF5Writer();

        virtual void initHDF5Groups();
        virtual void flushOutput();

    protected:

        HDF5File &hdf5File;

    private:

};

#endif
