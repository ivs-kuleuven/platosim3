#include "ThermoElasticDriftFromRedNoise.h"



/**
 * \brief function to acces the protected constructor, if there is no drift instance yet
 */
DriftGenerator* ThermoElasticDriftFromRedNoise::Instance(ConfigurationParameters &configParams, double readoutTimeBeforeNextExposure)
{
    if(_instance == 0)
    {
        _instance = new ThermoElasticDriftFromRedNoise(configParams, readoutTimeBeforeNextExposure);
    }
    return _instance;
}




/**
 * \brief Constructor
 * 
 * \param configParams The configuration parameters from the input parameters file
 *
 * \param readoutTimeBeforeNextExposure Duration of the readout that takes place before the next exposure can start
 */

ThermoElasticDriftFromRedNoise::ThermoElasticDriftFromRedNoise(ConfigurationParameters &configParams, double readoutTimeBeforeNextExposure)
: lastYaw(0.0), lastPitch(0.0), lastRoll(0.0), internalTime(0.0)
{
    // Set the configuration parameters

    configure(configParams, readoutTimeBeforeNextExposure);

    // Seed the random generator. The seed should have been set by configure().
    // Initialise the standard normal distribution with mu=0, and sigma=1.0.

    driftNoiseGenerator.seed(driftNoiseSeed);
    normalDistribution = normal_distribution<double>(0.0, 1.0);
}










/**
 * \brief Destructor
 */

ThermoElasticDriftFromRedNoise::~ThermoElasticDriftFromRedNoise()
{

}







/**
 * \brief Configure this object using the parameters from the input parameters file
 * 
 * \param configParams  The configuration parameters
 *
 * \param readoutTimeBeforeNextExposure Duration of the readout that takes place before the next exposure can start
 */

void ThermoElasticDriftFromRedNoise::configure(ConfigurationParameters &configParams, double readoutTimeBeforeNextExposure)
{
    // Note that the inputfile lists the drift RMS values in [arcsec]

    yawRMS          = deg2rad(configParams.getDouble("Telescope/DriftYawRms") / 3600.);     // [rad]         
    pitchRMS        = deg2rad(configParams.getDouble("Telescope/DriftPitchRms") / 3600.);   // [rad]         
    rollRMS         = deg2rad(configParams.getDouble("Telescope/DriftRollRms") / 3600.);    // [rad]         
    driftTimeScale  = configParams.getDouble("Telescope/DriftTimeScale");                   // [s]    
    driftNoiseSeed  = configParams.getLong("RandomSeeds/DriftSeed");

    // We determine the drift time interval as a fraction of the drift time scale.
    // so that the changes in (yaw, pitch, roll) can still be reliably tracked.

    driftTimeInterval = driftTimeScale / 20.0;

    // Set the internal time to the time of the first exposure of the series
    // The alternative, setting internalTime to 0., has the disadvantage that the jittering
    // needs to fast-forward a lot when beginExposureNr is very large, which is slow.
    
    int beginExposureNr = configParams.getInteger("ObservingParameters/BeginExposureNr");
    double exposureTime = configParams.getDouble("ObservingParameters/ExposureTime");

    internalTime = beginExposureNr * (exposureTime + readoutTimeBeforeNextExposure);
}









/**
 * \brief Get the next (yaw, pitch, roll) values using a Brownian motion model. These yaw, pitch,
 *        and roll values are with respect to the original pointing (at t=0) of the telescope, 
 *        NOT with respect to the last pointing.
 * 
 * \note Also during CCD readout, the telescope drifts, to the user needs to take this into
 *       account when passing 'timeInterval'.
 * 
 * \param time  The time point for which the (yaw, pitch, roll) is requested [s]
 * 
 * \return (newYaw, newPitch, newRoll)  [rad]
 */

tuple<double, double, double> ThermoElasticDriftFromRedNoise::getNextYawPitchRoll(double time)
{
    // If the time interval is zero, return the last computed values. 

    if (time == internalTime)
    {
        make_tuple(lastYaw, lastPitch, lastRoll);
    }

    // If the time interval is negative, complain

    if (time < internalTime)
    {
         string msg = "ThermoElasticDriftFromRedNoise: requested time < internal time: "
                      + to_string(time) + " < " + to_string(internalTime) + "Cannot rewind.";
         throw IllegalArgumentException(msg);
    }

    // Use bind() to get a shorter normal01() function to generate random numbers instead of 
    // the cumbersome normalDistribition(driftNoiseGenerator). Note: the std::ref() is needed, 
    // otherwise a copy is passed and the generator would always return the same number.

    auto normal01 = std::bind(normalDistribution, ref(driftNoiseGenerator));

    // The time step with which the yaw, pitch and roll will be iteratively updated.
    // Normally this time step is the driftTimeInterval, but if the user-given timeInterval
    // is actually smaller than, take the latter.

    double timeStep = min(time - internalTime, driftTimeInterval);

    // Initialise the (yaw, pitch, roll) values with the last computed ones

    double newYaw = lastYaw;
    double newPitch = lastPitch;
    double newRoll = lastRoll;

    // Move from the internal time to the user-given time in steps of 'timeStep', 
    // each time updating the yaw, pitch, and roll.

    int n = 0;
    while (internalTime + n * timeStep < time)
    {
        newYaw   = exp(-timeStep/driftTimeScale) * newYaw   + yawRMS   * sqrt(timeStep/driftTimeScale) * normal01();
        newPitch = exp(-timeStep/driftTimeScale) * newPitch + pitchRMS * sqrt(timeStep/driftTimeScale) * normal01();
        newRoll  = exp(-timeStep/driftTimeScale) * newRoll  + rollRMS  * sqrt(timeStep/driftTimeScale) * normal01();
    
        n++;
    }

    // In case that the user-given timeInterval cannot be covered with an integral number
    // of 'timeSteps', there is a small time interval left which still needs to be covered

    timeStep = time - internalTime - (n-1) * timeStep;

    newYaw   = exp(-timeStep/driftTimeScale) * newYaw   + yawRMS   * sqrt(timeStep/driftTimeScale) * normal01();
    newPitch = exp(-timeStep/driftTimeScale) * newPitch + pitchRMS * sqrt(timeStep/driftTimeScale) * normal01();
    newRoll  = exp(-timeStep/driftTimeScale) * newRoll  + rollRMS  * sqrt(timeStep/driftTimeScale) * normal01();


    // Save the (yaw, pitch, roll) values for the next request

    lastYaw = newYaw;
    lastPitch = newPitch;
    lastRoll = newRoll;

    // Update the internal clock

    internalTime = time;
    
    // That's it!

    return make_tuple(newYaw, newPitch, newRoll);
}








/**
 * \brief Return the heartbeat interval of this Red Noise drift generator
 * 
 * \details The heartbeat interval is the drift time interval which is set to a fraction of 
 *          the drift time scale. so that the changes in (yaw, pitch, roll) can still be 
 *          reliably tracked.
 *          
 * \return  heartbeatInterval [s]
 */

double ThermoElasticDriftFromRedNoise::getHeartbeatInterval()
{
    return driftTimeInterval;
}

