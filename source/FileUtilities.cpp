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

/**
 * @brief Return the last line in the file.
 * 
 * @param filename: Filename.
 * 
 * @return Last line in the file.
 */
string FileUtilities::getLastLine(const string &filename)
{
    Log.info("FileUtilities: Opening jitter file " + filename);

    ifstream fs(filename);
    string lastLine;

    if (fs.is_open())
    {
        fs.seekg(-1, std::ios_base::end);
        if (fs.peek() == '\n')
        {
            // Start searching for \n occurrences
            fs.seekg(-1, std::ios_base::cur);
            for (int i = fs.tellg(); i > 0; i--)
            {
                if (fs.peek() == '\n')
                {
                    // Found
                    fs.get();
                    break;
                }
                // Move one character back
                fs.seekg(i, std::ios_base::beg);
            }
        }
        getline(fs, lastLine);
    }
    else
    {
        string msg = "FileUtilities: Cannot open file " + filename;

        Log.error(msg);
        throw FileException(msg);
    }

    fs.close();

    return lastLine;
}

/**
 * Extract the last timepoint in the file.  Assumed is that the file consists of columns
 * and the first column is the time.
 * 
 * @param filename: Filename.
 * 
 * @return Last timepoint in the file.
 */
double FileUtilities::getLastTimePoint(const string &filename)
{
    string lastLine = getLastLine(filename);

    istringstream buffer(lastLine);
    vector<double> numbers((istream_iterator<double>(buffer)), istream_iterator<double>());
    
    return numbers[0];
}