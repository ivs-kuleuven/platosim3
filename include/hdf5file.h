
#ifndef HDF5FILE_H
#define HDF5FILE_H

#include <string>
#include "H5Cpp.h"

using namespace std;


class HDF5File
{
    public:

        HDF5File();
        HDF5File(string filename, bool overwrite=false);
        ~HDF5File();

        void open(string filename, bool overwrite=false);
        void close();
        
        void writeAttribute(string groupName, string attributeName, string attributeValue);
        void writeAttribute(string groupName, string attributeName, int attributeValue);
        void writeAttribute(string groupName, string attributeName, long attributeValue);
        void writeAttribute(string groupName, string attributeName, double attributeValue);

        void writeArray(string groupName, string arrayName, int *array);
        void writeArray(string groupName, string arrayName, float *array);
        void writeArray(string groupName, string arrayName, double *array);
        void writeArray(string groupName, string arrayName, double **array);
        
        void readArray(string groupName, string arrayName, float **array);


    protected:

        bool fileExists(string filename);

    private:

        H5::H5File *file;
        bool fileIsOpen;

};



#endif
