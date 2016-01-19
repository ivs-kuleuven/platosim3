
#ifndef HDF5FILE_H
#define HDF5FILE_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "H5Cpp.h"
#include "logger.h"
#include "memory.h"
#include "armadillo"

using namespace std;



class HDF5File
{
    public:

        HDF5File();
        HDF5File(string filename, bool overwrite=false);
        ~HDF5File();

        void open(string filename, bool overwrite=false);
        void close();
        
        void createGroup(string groupName);

        void writeAttribute(string groupName, string attributeName, string attributeValue);
        void writeAttribute(string groupName, string attributeName, int attributeValue);
        void writeAttribute(string groupName, string attributeName, long attributeValue);
        void writeAttribute(string groupName, string attributeName, double attributeValue);

        void writeArray(string groupName, string arrayName, int*    array, int size);
        void writeArray(string groupName, string arrayName, float*  array, int size);
        void writeArray(string groupName, string arrayName, double* array, int size);
        void writeArray(string groupName, string arrayName, arma::Mat<double>& A);


        void readArray(string groupName, string arrayName, arma::Mat<double>& A);


    protected:

        bool fileExists(string filename);


    private:

        H5::H5File *file;
        bool fileIsOpen;

};



#endif
