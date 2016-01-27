#include <iostream>

#include "Logger.h"
#include "Exceptions.h"

using std::string;

const char * IllegalArgumentException::what() const throw() {
    return string("IllegalArgumentException: " + message).c_str();
}

IllegalArgumentException::IllegalArgumentException(const char * msg) {
    message = msg;
    Log.error(what());
}

IllegalArgumentException::IllegalArgumentException(std::string msg) {
    message = msg;
    Log.error(what());
}

IllegalArgumentException::~IllegalArgumentException() throw() {}

