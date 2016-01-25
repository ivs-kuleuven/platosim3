#include <iostream>

#include "Exceptions.h"

using std::string;

const char * IllegalArgumentException::what() const throw() {
    return string("IllegalArgumentException: " + message).c_str();
}

IllegalArgumentException::IllegalArgumentException(const char * msg) {
    message = msg;
}

IllegalArgumentException::IllegalArgumentException(std::string msg) {
    message = msg;
}

IllegalArgumentException::~IllegalArgumentException() throw() {}

