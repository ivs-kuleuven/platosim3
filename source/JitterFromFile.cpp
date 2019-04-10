
#include "JitterFromFile.h"



/**
 * \brief Constructor
 * 
 * \param configParams The configuration parameters from the input parameters file
 *
 * \param readoutTimeBeforeNextExposure Duration of the readout that takes place before the next exposure can start
 */

JitterFromFile::JitterFromFile(ConfigurationParameters &configParams, double readoutTimeBeforeNextExposure)
{
    // Set the configuration parameters

    configure(configParams, readoutTimeBeforeNextExposure);

    // Open the jitter file, and read time yaw, pitch, roll time series.
    // The time is assumed to be in [s], pitch, yaw, and roll in [arcsec].
    // The path of the jitter file should have been set in configure().

    Log.info("JitterFromFile: Opening jitter file " + pathToJitterFile);
    Log.info("JitterFromFile: Reading jitter steps from t=" + to_string(beginTime) + " to t=" + to_string(endTime));

    ifstream jitterFile(pathToJitterFile);
    if (jitterFile.is_open())
    {
        double previousTime, previousYaw, previousPitch, previousRoll;
        string line;
        while (getline(jitterFile, line))
        {
            istringstream buffer(line);
            vector<double> numbers((istream_iterator<double>(buffer)), istream_iterator<double>());
            double time = numbers[0];
            
            // Only read the part of the jitter file that is relevant for this simulation.
            // This saves a lot of time when the file is large but the simulation is short.
            // If the time points in the jitter file do not exactly coincide with the time range
            // then include one point before, and one point after. This is needed for the 
            // interpolation afterwards.
            // E.g. If the time range is [10, 20], and the jitter time points are 
            //           ..., 9.5, 10.5, ..., 19.5, 20.5
            //      then the points 9.5 until 20.5 should be read. 
            //      If the time range is [10, 20], and the jitter time points are
            //           ..., 10, 11, 12, ..., 20, 21, ...
            //      then the points 10 to 21 are read.
            
            if (time < beginTime)
            {
                // If the *next* line is within the simulation time range, then we should also keep
                // the current one.
                
                previousTime  = time;
                previousYaw   = deg2rad(numbers[1]/3600.);
                previousPitch = deg2rad(numbers[2]/3600.);
                previousRoll  = deg2rad(numbers[3]/3600.);
                continue;
            }
            
            // Check if we are beyond the simulation time range. If so, stop reading the jitter file.
            // By checking the last element of 'timeFromFile' rather than 'time', we ensure that we always
            // read one point too many, as we intend.
            
            if (timeFromFile.size() != 0)
            {
                if (timeFromFile.back() > endTime) break;
            }

            // If we arrive here, we are within the simulation's time range. 
            // If it's our first point within the proper time range, also persist the previous point.
            // Note: if the very first line contains the time == beginTime, then previousTime is not
            //       defined. Hence below: time > beginTime, rather than time >= beginTime. 
            
            if ((timeFromFile.size() == 0) & (time > beginTime)) 
            {
                timeFromFile.push_back(previousTime);   // [s]  
                yaw.push_back(previousYaw);             // [arcsec] -> [rad]
                pitch.push_back(previousPitch);         // [arcsec] -> [rad]
                roll.push_back(previousRoll);           // [arcsec] -> [rad]
            }

            // Persist the current jitter steps in vectors.
            
            timeFromFile.push_back(time);                    // [s]  
            yaw.push_back(deg2rad(numbers[1]/3600.));        // [arcsec] -> [rad]
            pitch.push_back(deg2rad(numbers[2]/3600.));      // [arcsec] -> [rad]
            roll.push_back(deg2rad(numbers[3]/3600.));       // [arcsec] -> [rad]
        }

        jitterFile.close();

        // If there are less than two points, we can't even derive a heartbeat interval

        if (timeFromFile.size() < 2)
        {
            string msg = "JitterFromFile: Jitter file contains less than 2 time points in relevant time range";
            Log.error(msg);
            throw FileException(msg);
        }
        else
        {        
            Log.info("JitterFromFile: found " + to_string(timeFromFile.size()) + " time points in jitter input file in relevant time range");
        }
    }
    else
    {
        string msg = "JitterFromFile: Cannot open jitter file " + pathToJitterFile;
        Log.error(msg);
        throw FileException(msg);
    }

    // We start with the first time

    timeIndex = 0;
    internalTime = timeFromFile[0];

    // Log the heartbeat interval of this jitter generator
    
    double heartBeatInterval = getHeartbeatInterval();
    Log.info("JitterFromFile: heartbeat interval: " + to_string(heartBeatInterval));

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

void JitterFromFile::configure(ConfigurationParameters &configParams, double readoutTimeBeforeNextExposure)
{
    pathToJitterFile = configParams.getAbsoluteFilename("Platform/JitterFileName");
    int numExposures      = configParams.getInteger("ObservingParameters/NumExposures");
    int beginExposureNr   = configParams.getInteger("ObservingParameters/BeginExposureNr");
    double exposureTime   = configParams.getDouble("ObservingParameters/ExposureTime");

    //  Determine from when to when the simulation runs. Only for this time interval
    //  we need to read the jitter file into memory. This saves time when the jitter
    //  file is large but the simulation is short.
    
    beginTime = beginExposureNr * (exposureTime + readoutTimeBeforeNextExposure);
    endTime   = (beginExposureNr + numExposures) * (exposureTime + readoutTimeBeforeNextExposure);
}










/**
 * \brief Get the next (yaw, pitch, roll) values from the pre-computed jitter series.
 * 
 * \details A linear interpolation between time points will be done if the timeInterval 
 *          does not coincide with the time step of the jitter file.
 * 
 * \param time  The time point for which the (yaw, pitch, roll) is requested [s]
 * 
 * \return (newYaw, newPitch, newRoll)   [rad]
 */

tuple<double, double, double> JitterFromFile::getNextYawPitchRoll(double time)
{
    // If the time interval is zero then no interpolation is needed, just return 
    // the (yaw, pitch, roll) at the current index. Don't advance the timeIndex, 
    // because more timeInterval==0.0 may turn up.

    if (time == internalTime)
    {
        return make_tuple(yaw[timeIndex], pitch[timeIndex], roll[timeIndex]);
    }

    // If timeInterval is larger than zero, Advance the pointer 'timeIndex' in our precomputed 
    // jitter series such that we have
    //      time[index] <= internalTime + timeInterval < time[index+1]

    while (timeFromFile[timeIndex] < time)
    {
        timeIndex++;
        if (timeIndex >= timeFromFile.size())
        {
            Log.error("JitterFromFile: jitter file not large enough");
            exit(1);
        }
    }

    // Do a linear interpolation

    timeIndex--;
    const double weight1 = (time - timeFromFile[timeIndex]) / (timeFromFile[timeIndex+1] - timeFromFile[timeIndex]);
    const double weight2 = (timeFromFile[timeIndex+1] - time) / (timeFromFile[timeIndex+1] - timeFromFile[timeIndex]);
    const double newYaw   = yaw[timeIndex]   * weight2 + yaw[timeIndex+1]   * weight1;
    const double newPitch = pitch[timeIndex] * weight2 + pitch[timeIndex+1] * weight1;
    const double newRoll  = roll[timeIndex]  * weight2 + roll[timeIndex+1]  * weight1;
    
    // Update the internal time

    internalTime = time;

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
    return (timeFromFile[1] - timeFromFile[0]);
}

