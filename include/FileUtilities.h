#ifndef FILE_UTILITIES_H
#define FILE_UTILITIES_H

#include <string>
#include <fstream>

using namespace std;

class FileUtilities
{
public:
    static bool fileExists(const string &filename);
    static bool isRelative(const string &filename);
    static void remove(const string &filename);

};



#endif /* FileUtilities_hpp */
