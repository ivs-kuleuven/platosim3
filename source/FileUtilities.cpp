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




/**
 * PURPOSE: Check if a path is relative or absolute. An absolute path starts 
 *          with a '/' character, otherwise the path is considered relative.
 *
 * INPUTS:  path a path name
 *
 * OUTPUTS: true if path is relative, false otherwise
 */
bool FileUtilities::isRelative(const string &path) {
    if (path[0] == '/')
    {
        return false;
    }
    else
    {
        return true;
    }
}
