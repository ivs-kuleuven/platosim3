#include <iostream>

#include "Logger.h"
#include "Exceptions.h"

using std::string;

const char * FileException::what() const throw() {
    return string("FileException: " + message).c_str();
}

FileException::FileException(const char * msg) {
    message = msg;
    Log.error(what());
}

FileException::FileException(std::string msg) {
    message = msg;
    Log.error(what());
}

FileException::~FileException() throw() {}

