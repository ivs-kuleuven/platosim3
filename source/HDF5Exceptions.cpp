#include <iostream>

#include "Logger.h"
#include "HDF5Exceptions.h"

using std::string;





const char * H5FileException::what() const throw() {
    return string("H5FileException: " + message).c_str();
}

H5FileException::H5FileException(const char * msg) {
    message = msg;
    Log.error(what());
}

H5FileException::H5FileException(std::string msg) {
    message = msg;
    Log.error(what());
}

H5FileException::~H5FileException() throw() {}








const char * H5GroupException::what() const throw() {
    return string("H5GroupException: " + message).c_str();
}

H5GroupException::H5GroupException(const char * msg) {
    message = msg;
    Log.error(what());
}

H5GroupException::H5GroupException(std::string msg) {
    message = msg;
    Log.error(what());
}

H5GroupException::~H5GroupException() throw() {}








const char * H5DatasetException::what() const throw() {
    return string("H5DatasetException: " + message).c_str();
}

H5DatasetException::H5DatasetException(const char * msg) {
    message = msg;
    Log.error(what());
}

H5DatasetException::H5DatasetException(std::string msg) {
    message = msg;
    Log.error(what());
}

H5DatasetException::~H5DatasetException() throw() {}








const char * H5AttributeException::what() const throw() {
    return string("H5AttributeException: " + message).c_str();
}

H5AttributeException::H5AttributeException(const char * msg) {
    message = msg;
    Log.error(what());
}

H5AttributeException::H5AttributeException(std::string msg) {
    message = msg;
    Log.error(what());
}

H5AttributeException::~H5AttributeException() throw() {}

