//
//  Exceptions.h
//  New PLATO Simulator
//
//  Created by Rik Huygen on 27/11/15.
//  Copyright © 2015 KU Leuven. All rights reserved.
//

#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <string>
#include <exception>

class IOException : public std::exception {
    std::string message;
public:
    IOException(const char * msg);
    IOException(const std::string msg);
    virtual ~IOException() throw();
    const char * what() const throw();
};



class UnsupportedException : public std::exception {
    std::string message;
public:
    UnsupportedException(const char * msg);
    UnsupportedException(const std::string msg);
    virtual ~UnsupportedException() throw();
    const char * what() const throw();
};




#endif /* EXCEPTIONS_H */
