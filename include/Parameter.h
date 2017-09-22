
#ifndef PARAMETER_H
#define PARAMETER_H



#include <fstream>
#include <string>
#include <utility>
#include <sstream>
#include <iostream>
#include "Exceptions.h"


using namespace std;



template <typename T>
class Parameter
{
    public:
        
        Parameter(T fixedValue);
        Parameter(string filePath, T unitConversionFactor);
        ~Parameter();

        const T operator()(const double time);
        const T operator()();
        void updateValue(double time);
        bool isTimeDependent();

    protected:

        void openFile(string filePath);
        pair<double, T> readNextFromFile();

    private:

        bool isFixedToValue;         // True if the parameter is constant in time, false otherwise.
        T currentValue;              // The parameter value for the last updated time point 
        string filePath;             // Name of time series file, if the parameter is not constant in time
        ifstream inputFile;          // Time series file
        T unitConversionFactor;      // Unit conversion factor to multiply the values in the file with.

        double time1;                // Only two consecutive pairs (time1, value1) and (time2, value2)
        double time2;                //      are read and stored (rather than reading and storing the entire
        T value1;                    //      time series in memory at once), to save memory. As a consequence
        T value2;                    //      one cannot go 'back' in time and rewind the file reading.
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
Parameter<T>::Parameter(T fixedValue)
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
Parameter<T>::Parameter(string filePath, T unitConversionFactor)
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
Parameter<T>::~Parameter()
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

template <typename T>
void Parameter<T>::updateValue(double time)
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
    }


    const double weight1 = (time - time1) / (time2 - time1); 
    const double weight2 = (time2 - time) / (time2 - time1); 
    currentValue = value1 * weight2 + value2 * weight1;

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
const T Parameter<T>::operator()(const double time)
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
const T Parameter<T>::operator()()
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
bool Parameter<T>::isTimeDependent()
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
void Parameter<T>::openFile(string filePath)
{
    inputFile.open(filePath, ifstream::in);
    if (!inputFile.good() || !inputFile.is_open())
    {
        string msg = "Unable to open parameter time series from " + filePath;
        throw FileException(msg);
    }
}







/**
 * \brief Read the next line in the file and return a (time, value) pair.
 *
 * \param -
 *
 * \return A pair containing the time and the parameter value.
 */

template <typename T>
pair<double, T> Parameter<T>::readNextFromFile()
{
    pair<double, T> timeAndValue;

    string buffer;

    getline(inputFile, buffer);
    stringstream reader(buffer);
    reader >> timeAndValue.first;
    reader >> timeAndValue.second;

    // The values in the file might be SI units, while the code may need a different unit.
    // Hence we allow the user the specify a unit conversion factor. 

    timeAndValue.second *= unitConversionFactor;       

    // That's it!

    return timeAndValue;
}



#endif
