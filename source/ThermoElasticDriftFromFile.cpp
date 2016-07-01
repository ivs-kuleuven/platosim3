
#include "ThermoElasticDriftFromFile.h"



/**
 * \brief Constructor
 * 
 * \param configParams The configuration parameters from the input parameters file
 */

ThermoElasticDriftFromFile::ThermoElasticDriftFromFile(ConfigurationParameters &configParams)
{
    // Set the configuration parameters

    configure(configParams);

    // Open the thermo-elastic drift file, and read time yaw, pitch, roll time series.
    // The time is assumed to be in [s], pitch, yaw, and roll in [arcsec].
    // The path of the thermo-elastic drift file should have been set in configure().

    Log.info("ThermoElasticDriftFromFile: Opening input file " + pathToDriftFile + " to read");

    ifstream driftFile(pathToDriftFile);
    if (driftFile.is_open())
    {
        string temp;
        while (getline(driftFile, temp))
        {
            istringstream buffer(temp);
            vector<double> numbers((istream_iterator<double>(buffer)), istream_iterator<double>());
            time.push_back(numbers[0]);                      // [s]  
            yaw.push_back(deg2rad(numbers[1]/3600.));        // [arcsec] -> [rad]
            pitch.push_back(deg2rad(numbers[2]/3600.));      // [arcsec] -> [rad]
            roll.push_back(deg2rad(numbers[3]/3600.));       // [arcsec] -> [rad]
        }

        driftFile.close();

        // If there are less than two points, we can't even derive a heartbeat interval

        if (time.size() < 2)
        {
            Log.error("ThermoElasticDriftFromFile: Jitter file contains less than 2 time points");
            // FIXME: exit ???
        }
        else
        {        
            Log.info("ThermoElasticDriftFromFile: found " + to_string(time.size()) + " time points in jitter input file");
        }
    }
    else
    {
        Log.error("ThermoElasticDriftFromFile: Cannot open jitter file " + pathToDriftFile);
        exit(1);
    }

    // We start with the first time point of the file

    timeIndex = 0;
    internalTime = time[0];

}










/**
 * \brief Destructor
 */

ThermoElasticDriftFromFile::~ThermoElasticDriftFromFile()
{

}







/**
 * \brief Configure this object using the parameters from the input parameters file
 * 
 * \param configParams  The configuration parameters
 */

void ThermoElasticDriftFromFile::configure(ConfigurationParameters &configParams)
{
    pathToDriftFile = configParams.getAbsoluteFilename("Telescope/DriftFileName");
}










/**
 * \brief Get the next (yaw, pitch, roll) values from the pre-computed drift series.
 * 
 * \details A linear interpolation between time points will be done if the timeInterval 
 *          does not coincide with the time step of the drift file.
 * 
 * \param timeInterval[in]  Time interval that has passed since the last getNextYawPitchRoll() request. [s]
 * 
 * \return (newYaw, newPitch, newRoll)   [rad]
 */

tuple<double, double, double> ThermoElasticDriftFromFile::getNextYawPitchRoll(double timeInterval)
{
    // If the time interval is zero then no interpolation is needed, just return 
    // the (yaw, pitch, roll) at the current index. Don't advance the timeIndex, 
    // because more timeInterval==0.0 may turn up.

    if (timeInterval == 0.0)
    {
        return make_tuple(yaw[timeIndex], pitch[timeIndex], roll[timeIndex]);
    }

    // If timeInterval is larger than zero, Advance the pointer 'timeIndex' in our precomputed 
    // jitter series such that we have
    //      time[index] <= internalTime + timeInterval < time[index+1]

    while (time[timeIndex] < internalTime + timeInterval)
    {
        timeIndex++;
        if (timeIndex >= time.size())
        {
            Log.error("ThermoElasticDriftFromFile: input file not large enough");
            exit(1);
        }
    }

    // Do a linear interpolation

    timeIndex--;
    const double weight1 = (internalTime + timeInterval - time[timeIndex]) / (time[timeIndex+1] - time[timeIndex]);
    const double weight2 = (time[timeIndex+1] - internalTime - timeInterval) / (time[timeIndex+1] - time[timeIndex]);
    const double newYaw   = yaw[timeIndex]   * weight2 + yaw[timeIndex+1]   * weight1;
    const double newPitch = pitch[timeIndex] * weight2 + pitch[timeIndex+1] * weight1;
    const double newRoll  = roll[timeIndex]  * weight2 + roll[timeIndex+1]  * weight1;

    // Update the internal time

    internalTime += timeInterval;

    // That's it!
    
    return make_tuple(newYaw, newPitch, newRoll);
}









/**
 * \brief Return the heartbeat interval of this thermo-elastic drift generator
 * 
 * \details It is assumed that the pre-computed drift time series are equidistant so that
 *          the time step over which the drift changes is simply time[1]-time[0].
 *          
 * \return  heartbeatInterval [s]
 */

double ThermoElasticDriftFromFile::getHeartbeatInterval()
{
    return (time[1] - time[0]);
}

