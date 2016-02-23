#ifndef STRING_UTILITIES_H
#define STRING_UTILITIES_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cstddef>

using namespace std;

class StringUtilities
{
public:
    static bool ends_with(string const &, string const &);
    static vector<string> split(string, char);
    static void print(vector <string> &);
    static string dtos(double value, bool scientific = false, int precision = 6);

};


#endif /* STRING_UTILITIES_H */
