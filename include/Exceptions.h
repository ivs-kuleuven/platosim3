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



class IllegalArgumentException : public std::exception {
    std::string message;
public:
    IllegalArgumentException(const char * msg);
    IllegalArgumentException(const std::string msg);
    virtual ~IllegalArgumentException() throw();
    const char * what() const throw();
};




#endif /* EXCEPTIONS_H */
