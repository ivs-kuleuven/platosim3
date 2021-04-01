#include "TemperatureFromFile.h"



/**
 * \brief Constructor.
 *
 * \param configParams Configuration parameters from the input parameters file.
 *
 * \param component Component for which the temperature is stored in a file (FEE or CCD).
 */

TemperatureFromFile::TemperatureFromFile(ConfigurationParameters &configParams, string component)
{
    // Set the configuration parameters

    configure(configParams, component);

    // Open the temperature time series.
    // The time is assumed to be in [s], temperature in [K].
    // The path of the temperature file should have been set in configure().

    Log.info("TemperatureFromFile: Opening input file " + pathToTemperatureFile + " to read temperature for the " + component);

    ifstream temperatureFile(pathToTemperatureFile);

    if (temperatureFile.is_open())
    {
        string temp;
        while (getline(temperatureFile, temp))
        {
            istringstream buffer(temp);
            vector<double> numbers((istream_iterator<double>(buffer)), istream_iterator<double>());
            timeFromFile.push_back(numbers[0]);       // [s]
            temperature.push_back(numbers[1]);        // [K]
        }

        temperatureFile.close();

        if (timeFromFile.size() < 2)
        {
            string msg = "TemperatureFromFile: Temperature file for the " + component + " contains less than 2 time points";
            Log.error(msg);
            throw FileException(msg);
        }
        else
        {
            Log.info("TemperatureFromFile: found " + to_string(timeFromFile.size()) + " time points in temperature input file for the " + component);
        }
    }
    else
    {
        string msg = "TemperatureFromFile: Cannot open temperature file " + pathToTemperatureFile;
        Log.error(msg);
        throw FileException(msg);
    }

    // We start with the first time point of the file

    timeIndex = 0;
    internalTime = timeFromFile[0];
    currentTemperature = temperature[0];
}










/**
 * \brief Destructor.
 */

TemperatureFromFile::~TemperatureFromFile()
{

}










/**
 * \brief Configure the temperature variations for the given component, based on the given configuration parameters.
 *
 * \param configParams Configuration parameters.
 *
 * \param component Component for which the temperature variations over time are stored in a file.
 */

void TemperatureFromFile::configure(ConfigurationParameters &configParams, string component)
{
    pathToTemperatureFile = configParams.getAbsoluteFilename(component + "/TemperatureFileName");
}










/**
 * \brief Get the next temperature value from the pre-defined temperature time series.
 *
 * \details A linear interpolation between time points will be done if the timeInterval
 *          does not coincide with the time step of the temperature file.
 *
 * \param time  The time point for which the temperature is requested [s]
 *
 * \return Temperature [K]
 */
double TemperatureFromFile::getNextTemperature(double time)
{
    // If the time interval is zero then no interpolation is needed, just return
    // the temperature at the current index. Don't advance the timeIndex,
    // because more timeInterval==0.0 may turn up.

    if (time == internalTime)
    {
        return currentTemperature;
    }

    // If timeInterval is larger than zero, Advance the pointer 'timeIndex' in our pre-defined
    // temperature time series such that we have
    //      time[index] <= internalTime + timeInterval < time[index+1]

    while (timeFromFile[timeIndex] < time)
    {
        timeIndex++;
        if (timeIndex >= timeFromFile.size())
        {
            Log.error("TemperatureFromFile: input file " + pathToTemperatureFile + " not large enough");
            exit(1);
        }
    }

    // Do a linear interpolation

    timeIndex--;
    const double weight1 = (time - timeFromFile[timeIndex]) / (timeFromFile[timeIndex+1] - timeFromFile[timeIndex]);
    const double weight2 = (timeFromFile[timeIndex+1] - time) / (timeFromFile[timeIndex+1] - timeFromFile[timeIndex]);

    currentTemperature = temperature[timeIndex] * weight2 + temperature[timeIndex+1] * weight1;

    // Update the internal time

    internalTime = time;

    // That's it!

    return currentTemperature;
}
