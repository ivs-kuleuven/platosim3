
#ifndef PARAMETER_H
#define PARAMETER_H



#include <fstream>
#include <string>
#include <utility>
#include <sstream>
#include <iostream>
#include <array>
#include "Exceptions.h"


using namespace std;



// A Parameter may be a scalar or an array of scalars.
// // The template parametes are:
//   T : the type of the parameter (usually double)
//   N : the number of scalars the Parameter consist of (usually only 1)
//
// For N=1 (e.g. in case of the focal length) we don't want to deal with 
// e.g. array<double, 1> but prefer to use directly 'double'. To achieve
// this we need to use type-traits and use template specialisation of N=1.


template<typename T, unsigned int N>
struct valueType
{ 
    typedef array<T,N> type; 
};


template<typename T>
struct valueType<T, 1>
{ 
    typedef T type; 
};




// --------------------------------------------------------
// The general case: Parameter<T,N>
//
// The template-specialised case of N=1 follows underneath
//---------------------------------------------------------


template <typename T, unsigned int N=1>
class Parameter
{
    public:
        
        Parameter(typename valueType<T,N>::type &fixedValue);
        Parameter(string filePath, T unitConversionFactor);
        ~Parameter();

        const typename valueType<T,N>::type operator()(const double time);
        const typename valueType<T,N>::type operator()();
        void updateValue(double time);
        bool isTimeDependent();

    protected:

        void openFile(string filePath);
        pair<double, typename valueType<T,N>::type> readNextFromFile();

    private:

        bool isFixedToValue;                         // True if the parameter is constant in time, false otherwise.
        typename valueType<T,N>::type currentValue;  // The parameter value for the last updated time point 
        string filePath;                             // Name of time series file, if the parameter is not constant in time
        ifstream inputFile;                          // Time series file
        T unitConversionFactor;                      // Unit conversion factor to multiply the values in the file with.

        double time1;                                // Only two consecutive pairs (time1, value1) and (time2, value2)
        double time2;                                //   are read and stored (rather than reading and storing the entire
        typename valueType<T,N>::type value1;        //   time series in memory at once), to save memory. As a consequence
        typename valueType<T,N>::type value2;        //   one cannot go 'back' in time and rewind the file reading.
};








/**
 * \brief Using this constructor, the parameter is assumed to 
 *        be constant in time, with the given value.
 *
 * \param fixedValue: the constant value of the parameter.
 *
 * \return -
 */

template <typename T, unsigned int N>
Parameter<T,N>::Parameter(typename valueType<T,N>::type &fixedValue)
    : isFixedToValue(true), currentValue(fixedValue)
{

}







/**
 * \brief Using this constructor, the parameter is assumed to be time dependent
 *        with the relevant time series stored in an ascii file with exactly two
 *        columns, seperated with a space: time and parameter value.
 *
 * \param filePath: the absolute path of the time series ascii file.
 *
 * \return -
 */

template <typename T, unsigned int N>
Parameter<T,N>::Parameter(string filePath, T unitConversionFactor)
    : isFixedToValue(false), filePath(filePath), unitConversionFactor(unitConversionFactor)
{
   openFile(filePath);

   // Read the first two lines of the input file, so that there is something to interpolate with
   
   tie(time1, value1) = readNextFromFile();
   tie(time2, value2) = readNextFromFile();

   if (!inputFile.good())
   {
        string msg = "Parameter::updateValue: cannot read first two (time,value) pairs from file " + filePath;
        throw FileException(msg);
   }
   else
   {
       currentValue = value1;
   }
}







/**
 * \brief Desctructor
 *
 * \param -
 *
 * \return -
 */

template <typename T, unsigned int N>
Parameter<T,N>::~Parameter()
{
    if (!isFixedToValue) inputFile.close(); 
}







/**
 * \brief Fast forward the reading of the time series in the input ascii file
 *        until you find two time values time1 and time2 that bracket the value 
 *        time in the argument: time1 <= time < time2. Then, linearly interpolate 
 *        to find the value of the parameter at the given time.
 *
 * \param time: time value. Unit: same as in the time series file.
 *
 * \return -
 */

template <typename T, unsigned int N>
void Parameter<T,N>::updateValue(double time)
{
    if (isFixedToValue) return;

    if (time < time1)
    {
        string msg = "Parameter::updateValue: cannot go back in time in file " + filePath;
        msg += " (time = " + to_string(time) + " < " + to_string(time1) + ")";
        throw FileException(msg);
    }

    while(time > time2)
    {
        time1 = time2;
        value1 = value2;
        tie(time2, value2) = readNextFromFile();

        if (!inputFile.good())
        {
            string msg = "Parameter::updateValue: not enough time points in file " + filePath;
            msg += " (time = " + to_string(time) + " goes beyond end of file)";
            throw FileException(msg);
        }

        if (time1 == time2)
        {
            string msg = "Parameter::updateValue: Time " + to_string(time1) + " occurs more than once in file " + filePath;
            throw FileException(msg);
        }

        if (time2 < time1)
        {
            string msg = "Parameter::updateValue: time values (col 0) are not ordered ascending in file " + filePath;
            throw FileException(msg);
        }

    }


    const double weight1 = (time - time1) / (time2 - time1); 
    const double weight2 = (time2 - time) / (time2 - time1); 
    for (unsigned int n = 0; n < N; ++n)
    {
        currentValue[n] = value1[n] * weight2 + value2[n] * weight1;
    }
}









/**
 * \brief Shorthand for:
 *        1) updateTime(time)
 *        2) return the current parameter value.
 *        If the parameter is constant in time, then the argument 'time' is not used.
 *
 * \param time: the value for which the parameter value should be returned.
 *
 * \return currentValue: the parameter value.
 */

template <typename T, unsigned int N>
const typename valueType<T,N>::type Parameter<T,N>::operator()(const double time)
{
    updateValue(time);
    return currentValue;
}







/**
 * \brief Return the parameter value of the last updated value of time.
 *        In case a file with a time series is used, the value is initialised with the value 
 *        of the first point. 
 *        
 * \param -
 *
 * \return currentValue: current value of the parameter
 */

template <typename T, unsigned int N>
const typename valueType<T,N>::type Parameter<T,N>::operator()()
{
    return currentValue;
}










/**
 * \brief Return false if the parameter was set to a fixed value with its constructor.
 *        True otherwise.
 *
 * \param -
 *
 * \return See above.
 */

template <typename T, unsigned int N>
bool Parameter<T,N>::isTimeDependent()
{
    return !isFixedToValue;
}








/**
 * \brief Open the ascii file with the time series of the parameter.
 *
 * \param filePath: the absolute path of the time series file
 *
 * \return -
 */

template <typename T, unsigned int N>
void Parameter<T,N>::openFile(string filePath)
{
    inputFile.open(filePath, ifstream::in);
    if (!inputFile.good() || !inputFile.is_open())
    {
        string msg = "Parameter: Unable to open parameter time series from " + filePath;
        throw FileException(msg);
    }
}







/**
 * \brief Read the next line in the file and return a (time, value) pair.
 *        where 'value' is an array.
 *
 * \param -
 *
 * \return A pair containing the time and the parameter value.
 */

template <typename T, unsigned int N>
pair<double, typename valueType<T,N>::type> Parameter<T,N>::readNextFromFile()
{
    double time;
    typename valueType<T,N>::type value;

    string buffer;

    getline(inputFile, buffer);
    stringstream reader(buffer);
    reader >> time;

    for (unsigned int n = 0; n < N; ++n)
    {
        reader >> value[n];
    }

    // The values in the file might be SI units, while the code may need a different unit.
    // Hence we allow the user the specify a unit conversion factor. 

    for (unsigned int n = 0; n < N; ++n)
    {
        value[n] *= unitConversionFactor;       
    }

    // That's it!

    return make_pair(time, value);
}












//-------------------------------------------------
// The template-specialised case of Parameter<T,1>
// ------------------------------------------------





template <typename T>
class Parameter<T,1>
{
    public:
        
        Parameter(typename valueType<T,1>::type &fixedValue);
        Parameter(string filePath, T unitConversionFactor);
        ~Parameter();

        const typename valueType<T,1>::type operator()(const double time);
        const typename valueType<T,1>::type operator()();
        void updateValue(double time);
        bool isTimeDependent();

    protected:

        void openFile(string filePath);
        pair<double, typename valueType<T,1>::type> readNextFromFile();

    private:

        bool isFixedToValue;                         // typename valueType<T,1>::typerue if the parameter is constant in time, false otherwise.
        typename valueType<T,1>::type currentValue;  // typename valueType<T,1>::typehe parameter value for the last updated time point 
        string filePath;                             // Name of time series file, if the parameter is not constant in time
        ifstream inputFile;                          // typename valueType<T,1>::typeime series file
        T unitConversionFactor;                      // Unit conversion factor to multiply the values in the file with.

        double time1;                                // Only two consecutive pairs (time1, value1) and (time2, value2)
        double time2;                                //   are read and stored (rather than reading and storing the entire
        typename valueType<T,1>::type value1;        //   time series in memory at once), to save memory. As a consequence
        typename valueType<T,1>::type value2;        //   one cannot go 'back' in time and rewind the file reading.
};










/**
 * \brief Using this constructor, the parameter is assumed to 
 *        be constant in time, with the given value.
 *
 * \param fixedValue: the constant value of the parameter.
 *
 * \return -
 */

template <typename T>
Parameter<T,1>::Parameter(typename valueType<T,1>::type &fixedValue)
    : isFixedToValue(true), currentValue(fixedValue)
{

}







/**
 * \brief Using this constructor, the parameter is assumed to be time dependent
 *        with the relevant time series stored in an ascii file with exactly two
 *        columns, seperated with a space: time and parameter value.
 *
 * \param filePath: the absolute path of the time series ascii file.
 *
 * \return -
 */

template <typename T>
Parameter<T,1>::Parameter(string filePath, T unitConversionFactor)
    : isFixedToValue(false), filePath(filePath), unitConversionFactor(unitConversionFactor)
{
   openFile(filePath);

   // Read the first two lines of the input file, so that there is something to interpolate with
   
   tie(time1, value1) = readNextFromFile();
   tie(time2, value2) = readNextFromFile();

   if (!inputFile.good())
   {
        string msg = "Parameter::updateValue: cannot read first two (time,value) pairs from file " + filePath;
        throw FileException(msg);
   }
   else
   {
       currentValue = value1;
   }
}







/**
 * \brief Desctructor
 *
 * \param -
 *
 * \return -
 */

template <typename T>
Parameter<T,1>::~Parameter()
{
    if (!isFixedToValue) inputFile.close(); 
}









/**
 * \brief Shorthand for:
 *        1) updateTime(time)
 *        2) return the current parameter value.
 *        If the parameter is constant in time, then the argument 'time' is not used.
 *
 * \param time: the value for which the parameter value should be returned.
 *
 * \return currentValue: the parameter value.
 */

template <typename T>
const typename valueType<T,1>::type Parameter<T,1>::operator()(const double time)
{
    updateValue(time);
    return currentValue;
}







/**
 * \brief Return the parameter value of the last updated value of time.
 *        In case a file with a time series is used, the value is initialised with the value 
 *        of the first point. 
 *        
 * \param -
 *
 * \return currentValue: current value of the parameter
 */

template <typename T>
const typename valueType<T,1>::type Parameter<T,1>::operator()()
{
    return currentValue;
}










/**
 * \brief Return false if the parameter was set to a fixed value with its constructor.
 *        True otherwise.
 *
 * \param -
 *
 * \return See above.
 */

template <typename T>
bool Parameter<T,1>::isTimeDependent()
{
    return !isFixedToValue;
}








/**
 * \brief Open the ascii file with the time series of the parameter.
 *
 * \param filePath: the absolute path of the time series file
 *
 * \return -
 */

template <typename T>
void Parameter<T,1>::openFile(string filePath)
{
    inputFile.open(filePath, ifstream::in);
    if (!inputFile.good() || !inputFile.is_open())
    {
        string msg = "Unable to open parameter time series from " + filePath;
        throw FileException(msg);
    }
}











/**
 * \brief Fast forward the reading of the time series in the input ascii file
 *        until you find two time values time1 and time2 that bracket the value 
 *        time in the argument: time1 <= time < time2. Then, linearly interpolate 
 *        to find the value of the parameter at the given time.
 *
 * \param time: time value. Unit: same as in the time series file.
 *
 * \return -
 */

template <typename T>
void Parameter<T,1>::updateValue(double time)
{
    if (isFixedToValue) return;

    if (time < time1)
    {
        string msg = "Parameter::updateValue: cannot go back in time in file " + filePath;
        msg += " (time = " + to_string(time) + " < " + to_string(time1) + ")";
        throw FileException(msg);
    }

    while(time > time2)
    {
        time1 = time2;
        value1 = value2;
        tie(time2, value2) = readNextFromFile();
        
        if (!inputFile.good())
        {
            string msg = "Parameter::updateValue: not enough time points in file " + filePath;
            msg += " (time = " + to_string(time) + " goes beyond end of file)";
            throw FileException(msg);
        }

        if (time1 == time2)
        {
            string msg = "Parameter::updateValue: Time " + to_string(time1) + " occurs more than once in file " + filePath;
            throw FileException(msg);
        }

        if (time2 < time1)
        {
            string msg = "Parameter::updateValue: time values (col 0) are not ordered ascending in file " + filePath;
            throw FileException(msg);
        }
    }


    const double weight1 = (time - time1) / (time2 - time1); 
    const double weight2 = (time2 - time) / (time2 - time1); 
    currentValue = value1 * weight2 + value2 * weight1;
}








/**
 * \brief Read the next line in the file and return a (time, value) pair.
 *        This is a template specialisation in case the Parameter type is
 *        a scalar (int, double, ...) rather than an array.
 *
 * \param -
 *
 * \return A pair containing the time and the parameter value.
 *
 */

template <typename T>
pair<double, typename valueType<T,1>::type> Parameter<T,1>::readNextFromFile()
{
    double time;
    T value;

    string buffer;

    getline(inputFile, buffer);
    stringstream reader(buffer);
    reader >> time;
    reader >> value;

    // The values in the file might be SI units, while the code may need a different unit.
    // Hence we allow the user the specify a unit conversion factor. 

    value *= unitConversionFactor;       

    // That's it!

    return make_pair(time, value);
}





#endif



