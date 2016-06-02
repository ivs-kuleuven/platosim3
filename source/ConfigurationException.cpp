#include <iostream>

#include "Logger.h"
#include "Exceptions.h"

using std::string;

const char * ConfigurationException::what() const throw()
{
    return message.c_str();
}


ConfigurationException::ConfigurationException(std::string msg)
{
    message = msg;
    Log.error(message);
}


ConfigurationException::~ConfigurationException() throw() {}

