
#ifndef HDF5FILE_H
#define HDF5FILE_H

#include "H5Cpp.h"

using namespace std;


class Hdf5File
{
    public:

        Hdf5File();
        Hdf5File(string filename, bool overwrite=false);
        ~Hdf5File();

        void open(string filename, bool overwrite=false);
        void close();
        
        void writeAttribute(double attribute, string groupName, string attributeName);
        void write1DIntArray(int *array, string groupName, string arrayName);
        void write1DFloatArray(float *array, string groupName, string arrayName);
        void write1DDoubleArray(double *array, string groupName, string arrayName);
        void write2DFloatArray(double **array, string groupName, string arrayName);
        void read2DFloatArray(float **array, string groupName, string arrayName);


    protected:

        bool fileExists(string filename);

    private:

        H5::H5File *file;
        bool fileIsOpen;

};