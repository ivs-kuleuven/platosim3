#include <iostream>

#include "Logger.h"
#include "Exceptions.h"

using std::string;

const char * IllegalArgumentException::what() const throw()
{
    return message.c_str();
}


IllegalArgumentException::IllegalArgumentException(std::string msg)
{
    message = msg;
    Log.error(message);
}


IllegalArgumentException::~IllegalArgumentException() throw() {}

