#include "JitterFromRedNoise.h"



/**
 * \brief Constructor
 * 
 * \param configParams The configuration parameters from the input parameters file
 */

JitterFromRedNoise::JitterFromRedNoise(ConfigurationParameters &configParams)
: lastYaw(0.0), lastPitch(0.0), lastRoll(0.0)
{
    // Set the configuration parameters

    configure(configParams);

    // Seed the random generator. The seed should have been set by configure().
    // Initialise the standard normal distribution with mu=0, and sigma=1.0.

    jitterNoiseGenerator.seed(jitterNoiseSeed);
    normalDistribution = normal_distribution<double>(0.0, 1.0);
}










/**
 * \brief Destructor
 */

JitterFromRedNoise::~JitterFromRedNoise()
{

}







/**
 * \brief Configure this object using the parameters from the input parameters file
 * 
 * \param configParams  The configuration parameters
 */

void JitterFromRedNoise::configure(ConfigurationParameters &configParams)
{
    // Note that the inputfile lists the jitter RMS values in [arcsec]

    yawRMS          = deg2rad(configParams.getDouble("Platform/JitterYawRms") / 3600.);          
    pitchRMS        = deg2rad(configParams.getDouble("Platform/JitterPitchRms") / 3600.);        
    rollRMS         = deg2rad(configParams.getDouble("Platform/JitterRollRms") / 3600.);         
    jitterTimeScale = configParams.getDouble("Platform/JitterTimeScale"); 
    jitterNoiseSeed = configParams.getLong("RandomSeeds/JitterSeed");

    // We determine the jitter time interval as a fraction of the jitter time scale.
    // so that the changes in (yaw, pitch, roll) can still be reliably tracked.

    jitterTimeInterval = jitterTimeScale / 20.0;
}









/**
 * \brief Get the next (yaw, pitch, roll) values using a Brownian motion model. These yaw, pitch,
 *        and roll values are with respect to the original pointing (at t=0), NOT with respect to 
 *        the last pointing.
 * 
 * \note Also during CCD readout, the spacecraft jitters, to the user needs to take this into
 *       account when passing 'timeInterval'.
 * 
 * \param timeInterval[in]  Time interval that has passed since the last getNextYawPitchRoll() request. [s]
 * 
 * \return (newYaw, newPitch, newRoll)  [rad]
 */

tuple<double, double, double> JitterFromRedNoise::getNextYawPitchRoll(double timeInterval)
{
    // If the time interval is zero, return the last computed values

    if (timeInterval == 0.0)
    {
        make_tuple(lastYaw, lastPitch, lastRoll);
    }

    // Use bind() to get a shorter normal01() function to generate random numbers instead of 
    // the cumbersome normalDistribition(jitterNoiseGenerator). Note: the std::ref() is needed, 
    // otherwise a copy is passed and the generator would always return the same number.

    auto normal01 = std::bind(normalDistribution, ref(jitterNoiseGenerator));

    // The time step with which the yaw, pitch and roll will be iteratively updated.
    // Normally this time step is the jitterTimeInterval, but if the user-given timeInterval
    // is actually smaller than, take the latter.

    double timeStep = min(timeInterval, jitterTimeInterval);

    // Initialise the (yaw, pitch, roll) values with the last computed ones

    double newYaw = lastYaw;
    double newPitch = lastPitch;
    double newRoll = lastRoll;

    // Move through the user-given timeInterval in steps of 'timeStep', 
    // each time updating the yaw, pitch, and roll.

    int n = 0;
    while (n * timeStep < timeInterval)
    {
        newYaw   = exp(-timeStep/jitterTimeScale) * newYaw   + yawRMS   * sqrt(timeStep/jitterTimeScale) * normal01();
        newPitch = exp(-timeStep/jitterTimeScale) * newPitch + pitchRMS * sqrt(timeStep/jitterTimeScale) * normal01();
        newRoll  = exp(-timeStep/jitterTimeScale) * newRoll  + rollRMS  * sqrt(timeStep/jitterTimeScale) * normal01();
    
        n++;
    }

    // In case that the user-given timeInterval cannot be covered with an integral number
    // of 'timeSteps', there is a small time interval left which still needs to be covered

    timeStep = timeInterval - (n-1) * timeStep;

    newYaw   = exp(-timeStep/jitterTimeScale) * newYaw   + yawRMS   * sqrt(timeStep/jitterTimeScale) * normal01();
    newPitch = exp(-timeStep/jitterTimeScale) * newPitch + pitchRMS * sqrt(timeStep/jitterTimeScale) * normal01();
    newRoll  = exp(-timeStep/jitterTimeScale) * newRoll  + rollRMS  * sqrt(timeStep/jitterTimeScale) * normal01();


    // Save the (yaw, pitch, roll) values for the next request

    lastYaw = newYaw;
    lastPitch = newPitch;
    lastRoll = newRoll;
    
    // That's it!

    return make_tuple(newYaw, newPitch, newRoll);
}








/**
 * \brief Return the heartbeat interval of this Red Noise jitter generator
 * 
 * \details The heartbeat interval is the jitter time interval which is set to a fraction of 
 *          the jitter time scale. so that the changes in (yaw, pitch, roll) can still be 
 *          reliably tracked.
 *          
 * \return  heartbeatInterval [s]
 */

double JitterFromRedNoise::getHeartbeatInterval()
{
    return jitterTimeInterval;
}

