#ifndef JSON_FORMAT_H_
#define JSON_FORMAT_H_

#include <string>

#include "ConfigurationFormat.h"
#include "json.hpp"

using json = nlohmann::json;

class JSonFormat : public ConfigurationFormat
{
public:
    JSonFormat(const std::string &);
    ~JSonFormat();

    int getInteger(const std::string &);
    double getDouble(const std::string &);
    std::string getString(const std::string &);

private:
    json j_input;

};


#endif /* JSON_FORMAT_H_ */
