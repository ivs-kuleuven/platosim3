
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

    // Check whether the required time span is covered by the jitter file

    double lastTimePoint = FileUtilities::getLastTimePoint(pathToDriftFile);
    double requiredTimeRange = configParams.getDouble("ObservingParameters/CycleTime") * (configParams.getInteger("ObservingParameters/NumExposures") + configParams.getInteger("ObservingParameters/BeginExposureNr"));

    if (lastTimePoint < requiredTimeRange)
    {
        string msg = "ThermoElasticDriftFromFile: required time span of " + to_string(requiredTimeRange) + "s not covered by drift file";
        Log.error(msg);
        throw FileException(msg);
    }

    // Open the thermo-elastic drift file, and read time yaw, pitch, roll time series.
    // The time is assumed to be in [s], pitch, yaw, and roll in [arcsec].
    // The path of the thermo-elastic drift file should have been set in configure().

    Log.info("ThermoElasticDriftFromFile: Opening input file " + pathToDriftFile + " to read");
    Log.info("ThermoElasticDriftFromFile: Reading drift steps from t=" + to_string(beginTime) + " to t=" + to_string(endTime));

    ifstream driftFile(pathToDriftFile);
    if (driftFile.is_open())
    {
        double previousTime, previousYaw, previousPitch, previousRoll;
        string line;
        while (getline(driftFile, line))
        {
            istringstream buffer(line);
            vector<double> numbers((istream_iterator<double>(buffer)), istream_iterator<double>());
            double time = numbers[0];

            // Only read the part of the drift file that is relevant for this simulation.
            // This saves a lot of time when the file is large but the simulation is short.
            // If the time points in the drift file do not exactly coincide with the time range
            // then include one point before, and one point after.
            // E.g. if the time range is [10, 20], and the drift time points are 
            //      ..., 9.5, 10.5, ..., 19.5, 20.5
            //      then the points 9.5 until 20.5 should be read.
            
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

            // Check if we are beyond the simulation time range. If so, stop reading the drift file.
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
			
			// Persist the drift steps in vectors.

            timeFromFile.push_back(time);                    // [s]  
            yaw.push_back(deg2rad(numbers[1]/3600.));        // [arcsec] -> [rad]
            pitch.push_back(deg2rad(numbers[2]/3600.));      // [arcsec] -> [rad]
            roll.push_back(deg2rad(numbers[3]/3600.));       // [arcsec] -> [rad]
        }

        driftFile.close();

        // If there are less than two points, we can't even derive a heartbeat interval

        if (timeFromFile.size() < 2)
        {
            string msg = "ThermoElasticDriftFromFile: Jitter file contains less than 2 time points in relevant time range";
            Log.error(msg);
            throw FileException(msg);
        }
        else
        {        
            Log.info("ThermoElasticDriftFromFile: found " + to_string(timeFromFile.size()) + " time points in drift input file in relevant time range");
        }
    }
    else
    {
        string msg = "ThermoElasticDriftFromFile: Cannot open drift file " + pathToDriftFile;
        Log.error(msg);
        throw FileException(msg);
    }

    // We start with the first time point of the file

    timeIndex = 0;
    internalTime = timeFromFile[0];

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
    int numExposures      = configParams.getInteger("ObservingParameters/NumExposures");
    int beginExposureNr   = configParams.getInteger("ObservingParameters/BeginExposureNr");
    double cycleTime   = configParams.getDouble("ObservingParameters/CycleTime");

    //  Determine from when to when the simulation runs. Only for this time interval
    //  we need to read the drift file into memory. This saves time when the drift
    //  file is large but the simulation is short.
    
    beginTime = beginExposureNr * cycleTime;
    endTime   = (beginExposureNr + numExposures) * cycleTime;
}










/**
 * \brief Get the next (yaw, pitch, roll) values from the pre-computed drift series.
 * 
 * \details A linear interpolation between time points will be done if the timeInterval 
 *          does not coincide with the time step of the drift file.
 * 
 * \param time  The time point for which the (yaw, pitch, roll) is requested [s]
 * 
 * \return (newYaw, newPitch, newRoll)   [rad]
 */

tuple<double, double, double> ThermoElasticDriftFromFile::getNextYawPitchRoll(double time)
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
            Log.error("ThermoElasticDriftFromFile: input file not large enough");
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
 * \brief Return the heartbeat interval of this thermo-elastic drift generator
 * 
 * \details It is assumed that the pre-computed drift time series are equidistant so that
 *          the time step over which the drift changes is simply time[1]-time[0].
 *          
 * \return  heartbeatInterval [s]
 */

double ThermoElasticDriftFromFile::getHeartbeatInterval()
{
    return (timeFromFile[1] - timeFromFile[0]);
}

