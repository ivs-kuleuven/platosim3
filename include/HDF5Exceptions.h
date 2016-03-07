#ifndef HDF5_EXCEPTIONS_H
#define HDF5_EXCEPTIONS_H

#include <string>
#include <exception>


using namespace std;



/**
 * \class H5FileException
 * \brief throw this exception for any problems with an HDF5 at file level
 */
class H5FileException : public exception 
{
    public:
        H5FileException(const string msg);
        virtual ~H5FileException() throw();
        virtual const char * what() const throw() override;
    private:
        string message;
};






/**
 * \class H5GroupException
 * \brief throw this exception for any problems with an HDF5 at group level
 */
class H5GroupException : public exception 
{
    public:
        H5GroupException(const string msg);
        virtual ~H5GroupException() throw();
        virtual const char * what() const throw() override;
    private:
        string message;
};






/**
 * \class H5DatasetException
 * \brief throw this exception for any problems with an HDF5 at dataset level
 */
class H5DatasetException : public exception 
{
    public:
        H5DatasetException(const string msg);
        virtual ~H5DatasetException() throw();
        virtual const char * what() const throw() override;
    private:
        string message;
};






/**
 * \class H5AttributeException
 * \brief throw this exception for any problems with an HDF5 at attribute level
 */
class H5AttributeException : public exception 
{
    public:
        H5AttributeException(const string msg);
        virtual ~H5AttributeException() throw();
        virtual const char * what() const throw() override;
    private:
        string message;
};


#endif /* HDF5_EXCEPTIONS_H */
