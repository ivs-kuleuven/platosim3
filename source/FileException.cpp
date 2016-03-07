#include <iostream>

#include "Logger.h"
#include "Exceptions.h"

using std::string;

const char * FileException::what() const throw()
{
    return message.c_str();
}


FileException::FileException(std::string msg)
{
    message = msg;
    Log.error(message);
}


FileException::~FileException() throw() {}

