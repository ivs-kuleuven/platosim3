#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <string>
#include <exception>



/**
 * \class IOException
 * \brief throw this exception for any input/output problems in your code
 */
class IOException : public std::exception {
    std::string message;
public:
    IOException(const char * msg);
    IOException(const std::string msg);
    virtual ~IOException() throw();
    const char * what() const throw();
};



/**
 * \class UnsupportedException
 * \brief throw this exception when some feature has not yet been supported in code
 */
class UnsupportedException : public std::exception {
    std::string message;
public:
    UnsupportedException(const char * msg);
    UnsupportedException(const std::string msg);
    virtual ~UnsupportedException() throw();
    const char * what() const throw();
};



/**
 * \class IllegalArgumentException
 * \brief throw this exception when an argument that was passed to a method or function is invalid
 */
class IllegalArgumentException : public std::exception {
    std::string message;
public:
    IllegalArgumentException(const char * msg);
    IllegalArgumentException(const std::string msg);
    virtual ~IllegalArgumentException() throw();
    const char * what() const throw();
};



/**
 * \class FileException
 * \brief throw this exception for any problems with files in your code
 */
class FileException : public std::exception {
    std::string message;
public:
    FileException(const char * msg);
    FileException(const std::string msg);
    virtual ~FileException() throw();
    const char * what() const throw();
};


#endif /* EXCEPTIONS_H */
