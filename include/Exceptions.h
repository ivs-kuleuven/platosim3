#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <string>
#include <exception>

using namespace std;



/**
 * \class IOException
 * \brief throw this exception for any input/output problems in your code
 */
class IOException : public exception 
{
    public:
        IOException(const string msg);
        virtual ~IOException() throw();
        virtual const char * what() const throw() override;
    private:
        string message;
};





/**
 * \class UnsupportedException
 * \brief throw this exception when some feature has not yet been supported in code
 */
class UnsupportedException : public exception 
{
    public:
        UnsupportedException(const string msg);
        virtual ~UnsupportedException() throw();
        virtual const char * what() const throw() override;
    private:
        string message;
};





/**
 * \class IllegalArgumentException
 * \brief throw this exception when an argument that was passed to a method or function is invalid
 */
class IllegalArgumentException : public exception 
{
    public:
        IllegalArgumentException(const string msg);
        virtual ~IllegalArgumentException() throw();
        virtual const char * what() const throw() override;
    private:
        string message;
};








/**
 * \class FileException
 * \brief throw this exception for any problems with files in your code
 */
class FileException : public exception 
{
    public:
        FileException(const string msg);
        virtual ~FileException() throw();
        virtual const char * what() const throw() override;
    private:
        string message;
};









/**
 * \class ConfigurationException
 * \brief throw this exception for any problems with the configuration, e.g. YAML
 */
class ConfigurationException : public exception 
{
    public:
        ConfigurationException(const string msg);
        virtual ~ConfigurationException() throw();
        virtual const char * what() const throw() override;
    private:
        string message;
};


#endif /* EXCEPTIONS_H */
