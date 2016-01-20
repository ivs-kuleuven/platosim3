#include <iostream>
#include "Exceptions.h"

using std::string;

const char * IOException::what() const throw() {
    return string("IOException: " + message).c_str();
}

IOException::IOException(const char * msg) {
    message = msg;
}

IOException::IOException(std::string msg) {
    message = msg;
}

IOException::~IOException() throw() {}

