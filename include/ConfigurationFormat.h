#ifndef CONFIGURATION_FORMAT_H_
#define CONFIGURATION_FORMAT_H_

#include <string>

class ConfigurationFormat
{
public:
    ConfigurationFormat();
    ~ConfigurationFormat();

    virtual int getInteger(const std::string &) = 0;
    virtual double getDouble(const std::string &) = 0;
    virtual std::string getString(const std::string &) = 0;

//    virtual int getIntegerArray(const std::string &);
//    virtual int getDoubleArray(const std::string &);
//    virtual int getStringArray(const std::string &);

};


#endif /* CONFIGURATION_FORMAT_H_ */