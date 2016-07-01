
#include "JitterFromFile.h"



/**
 * \brief Constructor
 * 
 * \param configParams The configuration parameters from the input parameters file
 */

JitterFromFile::JitterFromFile(ConfigurationParameters &configParams)
{
    // Set the configuration parameters

    configure(configParams);

    // Open the jitter file, and read time yaw, pitch, roll time series.
    // The time is assumed to be in [s], pitch, yaw, and roll in [arcsec].
    // The path of the jitter file should have been set in configure().

    Log.info("JitterFromFile: Opening jitter file " + pathToJitterFile + " to read");

    ifstream jitterFile(pathToJitterFile);
    if (jitterFile.is_open())
    {
        string temp;
        while (getline(jitterFile, temp))
        {
            istringstream buffer(temp);
            vector<double> numbers((istream_iterator<double>(buffer)), istream_iterator<double>());
            time.push_back(numbers[0]);                      // [s]  
            yaw.push_back(deg2rad(numbers[1]/3600.));        // [arcsec] -> [rad]
            pitch.push_back(deg2rad(numbers[2]/3600.));      // [arcsec] -> [rad]
            roll.push_back(deg2rad(numbers[3]/3600.));       // [arcsec] -> [rad]
        }

        jitterFile.close();

        // If there are less than two points, we can't even derive a heartbeat interval

        if (time.size() < 2)
        {
            Log.error("JitterFromFile: Jitter file contains less than 2 time points");
            // FIXME: exit ???
        }
        else
        {        
            Log.info("JitterFromFile: found " + to_string(time.size()) + " time points in jitter input file");
        }
    }
    else
    {
        Log.error("JitterFromFile: Cannot open jitter file " + pathToJitterFile);
        exit(1);
    }

    // We start with the first time

    timeIndex = 0;
    internalTime = time[0];

}










/**
 * \brief Destructor
 */

JitterFromFile::~JitterFromFile()
{

}







/**
 * \brief Configure this object using the parameters from the input parameters file
 * 
 * \param configParams  The configuration parameters
 */

void JitterFromFile::configure(ConfigurationParameters &configParams)
{
    pathToJitterFile = configParams.getAbsoluteFilename("Platform/JitterFileName");
}










/**
 * \brief Get the next (yaw, pitch, roll) values from the pre-computed jitter series.
 * 
 * \details A linear interpolation between time points will be done if the timeInterval 
 *          does not coincide with the time step of the jitter file.
 * 
 * \param timeInterval[in]  Time interval that has passed since the last getNextYawPitchRoll() request. [s]
 * 
 * \return (newYaw, newPitch, newRoll)   [rad]
 */

tuple<double, double, double> JitterFromFile::getNextYawPitchRoll(double timeInterval)
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
            Log.error("JitterFromFile: jitter file not large enough");
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
 * \brief Return the heartbeat interval of this Jitter generator
 * 
 * \details It is assumed that the pre-computed jitter time series are equidistant so that
 *          the time step over which the jitter changes is simply time[1]-time[0].
 *          
 * \return  heartbeatInterval [s]
 */

double JitterFromFile::getHeartbeatInterval()
{
    return (time[1] - time[0]);
}

