//
//  IOException.cpp
//  New PLATO Simulator
//
//  Created by Rik Huygen on 27/11/15.
//  Copyright © 2015 KU Leuven. All rights reserved.
//

#include <iostream>
#include "Exceptions.h"

using std::string;

const char * UnsupportedException::what() const throw() {
    return string("UnsupportedException: " + message).c_str();
}

UnsupportedException::UnsupportedException(const char * msg) {
    message = msg;
}

UnsupportedException::UnsupportedException(std::string msg) {
    message = msg;
}

UnsupportedException::~UnsupportedException() throw() {}

