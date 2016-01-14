//
//  FileUtilities.h
//  New PLATO Simulator
//
//  Created by Rik Huygen on 27/11/15.
//  Copyright © 2015 KU Leuven. All rights reserved.
//

#ifndef FILE_UTILITIES_H
#define FILE_UTILITIES_H

#include <string>
#include <sys/stat.h>

/**
 * PURPOSE: Check if a file exists
 * 
 * INPUTS:  filename
 * 
 * OUTPUTS: return true if and only if the file exists, false otherwise
 */
static bool fileExists(const std::string& file) {
    struct stat buf;
    return (stat(file.c_str(), &buf) == 0);
}

#endif /* FileUtilities_hpp */
