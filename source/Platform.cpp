
#include "Platform.h"



/**
 * \brief Constructor
 * 
 * \param configurationParameters  Configuration parameters as read from the (e.g. yaml) inputfile
 * \param hdf5File                 HDF5 file where to write any output
 * \param jitterGenerator          Generator of the yaw, pitch roll variations due to spacecraft jitter
 */

Platform::Platform(ConfigurationParameters configurationParameters, HDF5File &hdf5File, JitterGenerator &jitterGenerator)
: HDF5Writer(hdf5File), jitterGenerator(jitterGenerator)
{

}







/**
 * \brief Destructor
 */

Platform::~Platform()
{

}









void Platform::configure(ConfigurationParameters &configParams)
{
    // Currently empty. Jitter configuration is done through JitterGenerator.
}








/**
 * \brief Repoints the spacecraft to the given pointing coordinates
 * 
 * \param rightAscencsion  Equatorial right Ascension [see unit]  
 * \param declination      Equatorial declination [see unit] 
 * \param unit             Unit of the input coordinates. Either Angle::degrees or Angle::radians.
 */

void Platform::setPointingCoordinates(double rightAscencsion, double declination, Unit unit)
{
    currentRA = rightAscencsion / unit;
    currentDec = declination / unit;
}









/**
 * \brief Return the current pointing coordinates of the spacecraft
 * 
 * \return (alpha, delta)    RA & dec [rad] of the current pointing
 */

pair<double, double> Platform::getPointingCoordinates()
{
    return make_pair(currentRA, currentDec);
}









/**
 * \brief Update the pointing coordinates (varying due to jitter) until the given time point.
 * 
 * \param time [s]
 */

void Platform::updatePointingCoordinates(double time)
{
    internalTime = time;
    return;
}









/**
 * \brief Return the heartbeat interval of Platform. It's the largest interval with which 
 *        small Platform variations (e.g. due to Jitter) can be tracked.
 * 
 */

double Platform::getHeartbeatInterval()
{
    // Currently Platform does not change in any other way than just jittering, so delegate to the jitter generator.

    return jitterGenerator.getHeartbeatInterval();
}





