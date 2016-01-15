#include "FileUtilities.h"

/**
 * PURPOSE: Check if a file exists
 * 
 * INPUTS:  filename
 * 
 * OUTPUTS: return true if and only if the file exists, false otherwise
 */
bool FileUtilities::fileExists(const string &filename) {
    ifstream fin(filename);
    return fin.good();
}





bool FileUtilities::isRelative(const string &filename) {
    if (filename[0] == '/')
    {
        return false;
    }
    else
    {
        return true;
    }
}
