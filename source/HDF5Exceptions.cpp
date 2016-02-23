#include <iostream>

#include "Logger.h"
#include "HDF5Exceptions.h"

using namespace std;





const char * H5FileException::what() const throw() 
{
    return message.c_str();
}

H5FileException::H5FileException(string msg) 
{
    message = msg;
    Log.error(message);
}

H5FileException::~H5FileException() throw() {}








const char * H5GroupException::what() const throw() 
{
    return message.c_str();
}

H5GroupException::H5GroupException(string msg) 
{
    message = msg;
    Log.error(message);
}

H5GroupException::~H5GroupException() throw() {}








const char * H5DatasetException::what() const throw() 
{
    return message.c_str();
}

H5DatasetException::H5DatasetException(string msg) 
{
    message = msg;
    Log.error(message);
}

H5DatasetException::~H5DatasetException() throw() {}








const char * H5AttributeException::what() const throw() 
{
    return message.c_str();
}

H5AttributeException::H5AttributeException(string msg) 
{
    message = msg;
    Log.error(message);
}

H5AttributeException::~H5AttributeException() throw() {}



