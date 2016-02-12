
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

    // Open the jitter file, and read time yaw, pitch, roll time series
    // The path of the jitter file should have been set in configure().

    ifstream jitterFile(pathToJitterFile);
    if (jitterFile.is_open())
    {
        string temp;
        long n = 0;
        while (getline(jitterFile, temp))
        {
            istringstream buffer(temp);
            vector<double> numbers((istream_iterator<double>(buffer)), istream_iterator<double>());
            time[n]  = numbers[0];
            yaw[n]   = numbers[1];
            pitch[n] = numbers[2];
            roll[n]  = numbers[3];
            n++;
        }

        jitterFile.close();

        // If there are less than two points, we can't even derive a heartbeat interval

        if (time.size() < 2)
        {
            Log.error("JitterFromFile: Jitter file " + pathToJitterFile + " contains less than 2 time points");
        }
        else
        {        
            Log.info("JitterFromFile: found " + to_string(time.size()) + " time points in input file " + pathToJitterFile);
        }
    }
    else
    {
        Log.error("JitterFromFile: Cannot open jitter file " + pathToJitterFile);
        exit(1);
    }

    // We start with the first time

    timeIndex = 0;

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
    string rootPath       = configParams.getString("General/ProjectLocation");
    string jitterFileName = configParams.getString("Platform/JitterFileName");
    pathToJitterFile = rootPath + "/" + jitterFileName;
}









/**
 * \brief Get the next (yaw, pitch, roll) values from the pre-computed jitter series.
 * \details A linear interpolation between time points will be done if the timeInterval 
 *          does not coincide with the time step of the jitter file.
 * 
 * \param newYaw[out]       Will contain the next yaw value [?]
 * \param newPitch[out]     Will contain the next pitch value [?]
 * \param newRoll[out]      Will contain the next roll value [?]
 * \param timeInterval[in]  Time interval that has passed since the last getNextYawPitchRoll() request. [s]
 */

void JitterFromFile::getNextYawPitchRoll(double &newYaw, double &newPitch, double &newRoll, double timeInterval)
{
    // Advance the pointer 'timeIndex' in our precomputed jitter series such that we have
    //      time[index]-lastTime <= timeInterval < time[index+1]-lastTime 

    double lastTime = time[timeIndex];
    while (time[timeIndex] - lastTime < timeInterval)
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
    double weight1 = (lastTime + timeInterval - time[timeIndex]) / (time[timeIndex+1] - time[timeIndex]);
    double weight2 = (time[timeIndex+1] - lastTime - timeInterval) / (time[timeIndex+1] - time[timeIndex]);
    newYaw   = yaw[timeIndex]   * weight1 + yaw[timeIndex+1]   * weight2;
    newPitch = pitch[timeIndex] * weight1 + pitch[timeIndex+1] * weight2;
    newRoll  = roll[timeIndex]  * weight1 + roll[timeIndex+1]  * weight2;
    
    return;
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

