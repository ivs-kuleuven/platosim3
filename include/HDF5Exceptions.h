#ifndef HDF5_EXCEPTIONS_H
#define HDF5_EXCEPTIONS_H

#include <string>
#include <exception>


using namespace std;

/**
 * \class H5FileException
 * \brief throw this exception for any problems with an HDF5 at file level
 */
class H5FileException : public exception {
    string message;
public:
    H5FileException(const string msg);
    virtual ~H5FileException() throw();
    const char * what() const throw();
};






/**
 * \class H5GroupException
 * \brief throw this exception for any problems with an HDF5 at group level
 */
class H5GroupException : public exception {
    string message;
public:
    H5GroupException(const string msg);
    virtual ~H5GroupException() throw();
    const char * what() const throw();
};






/**
 * \class H5DatasetException
 * \brief throw this exception for any problems with an HDF5 at dataset level
 */
class H5DatasetException : public exception {
    string message;
public:
    H5DatasetException(const string msg);
    virtual ~H5DatasetException() throw();
    const char * what() const throw();
};






/**
 * \class H5AttributeException
 * \brief throw this exception for any problems with an HDF5 at attribute level
 */
class H5AttributeException : public exception {
    string message;
public:
    H5AttributeException(const string msg);
    virtual ~H5AttributeException() throw();
    const char * what() const throw();
};


#endif /* HDF5_EXCEPTIONS_H */
