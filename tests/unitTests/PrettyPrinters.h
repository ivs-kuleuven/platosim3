#ifndef PRETTY_PRINTERS_H
#define PRETTY_PRINTERS_H

#include <iostream>
#include <string>

#include "Logger.h"
#include "Exceptions.h"
#include "HDF5Exceptions.h"

// It's important that PrintTo() for each class is defined in the SAME
// namespace that defines the class.  C++'s look-up rules rely on that.




void PrintTo(const FileException& fe, ::std::ostream* os) {
    *os << fe.what();
}

void PrintTo(const IOException& ioe, ::std::ostream* os) {
    *os << ioe.what();
}

void PrintTo(const UnsupportedException& ue, ::std::ostream* os) {
    *os << ue.what();
}

void PrintTo(const IllegalArgumentException& iae, ::std::ostream* os) {
    *os << iae.what();
}

void PrintTo(const H5FileException& iae, ::std::ostream* os) {
    *os << iae.what();
}

void PrintTo(const H5GroupException& iae, ::std::ostream* os) {
    *os << iae.what();
}

void PrintTo(const H5DatasetException& iae, ::std::ostream* os) {
    *os << iae.what();
}

void PrintTo(const H5AttributeException& iae, ::std::ostream* os) {
    *os << iae.what();
}




#endif // PRETTY_PRINTERS_H