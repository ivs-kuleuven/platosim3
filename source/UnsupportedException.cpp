#include <iostream>

#include "Logger.h"
#include "Exceptions.h"

using std::string;

const char * UnsupportedException::what() const throw() {
    return string("UnsupportedException: " + message).c_str();
}

UnsupportedException::UnsupportedException(const char * msg) {
    message = msg;
    Log.error(what());
}

UnsupportedException::UnsupportedException(std::string msg) {
    message = msg;
    Log.error(what());
}

UnsupportedException::~UnsupportedException() throw() {}

