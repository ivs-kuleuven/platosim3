
#ifndef HDF5WRITER_H
#define HDF5WRITER_H

#include <limits>
#include <string>
#include "hdf5file.h"

using namespace std;



class Hdf5Writer
{
    
    public:

        Hdf5Writer(Hdf5File &hdf5File, string groupName);
        ~Hdf5Writer();

        virtual void flushOutput();

    protected:


    private:

        Hdf5File &hdf5File;
        string hdf5GroupName;
    
};

#endif
