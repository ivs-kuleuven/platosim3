#include <iostream>

#include "Logger.h"
#include "Exceptions.h"

using std::string;

const char * UnsupportedException::what() const throw()
{
    return message.c_str();
}


UnsupportedException::UnsupportedException(std::string msg)
{
    message = msg;
    Log.error(message);
}


UnsupportedException::~UnsupportedException() throw() {}

