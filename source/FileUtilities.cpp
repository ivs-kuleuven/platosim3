#include <cstdio>   // for remove file

#include "FileUtilities.h"
#include "Logger.h"

/**
 * @brief      Check if a file exists
 * 
 * @param[in]  filename
 * 
 * @returns    true if and only if the file exists, false otherwise
 */
bool FileUtilities::fileExists(const string &filename) {
    ifstream fin(filename);
    return fin.good();
}




/**
 * @brief      Check if a path is relative or absolute.
 * 
 * @details
 * 
 * An absolute path starts with a '/' character, otherwise the path is considered relative.
 *
 * @param[in]  path a path name
 *
 * @returns    true if path is relative, false otherwise
 * 
 * @todo       make this also work with Window paths, e.g. C:\\Documents
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



void FileUtilities::remove(const string &filename)
{
    if (fileExists(filename))
    {
        if (std::remove(filename.c_str()))
        {
            Log.warning("Couldn't remove the file: " + filename);
        }
    }
}