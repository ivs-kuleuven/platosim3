#include <iostream>

#include "Logger.h"
#include "Exceptions.h"

using std::string;

const char * IOException::what() const throw() 
{
    return message.c_str();
}


IOException::IOException(std::string msg)
{
    message = msg;
    Log.error(message);
}

IOException::~IOException() throw() {}

