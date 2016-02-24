#ifndef STRING_UTILITIES_H
#define STRING_UTILITIES_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cstddef>

using namespace std;



namespace StringUtilities
{
    bool ends_with(string const & value, string const & ending);
    vector<string> split(string myString, char delimiter);
    string dtos(double value, bool scientific = false, int precision = 6);
    void print( std::vector <std::string> & vector );
}



#endif /* STRING_UTILITIES_H */
