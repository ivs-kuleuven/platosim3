#include "Telescope.h"

/**
 * Constructor
 * 
 * \param configurationParameters: Configuration parameters for the telescope.
 * \param hdf5File                 Output HDF5 file.
 * \param Platform:                Platform on which the telescope is mounted
 * 
 */

Telescope::Telescope(ConfigurationParameters &configParams, HDF5File &hdf5File, Platform &platform)
: HDF5Writer(hdf5File), internalTime(0.0), platform(platform)
{
	// Retrieve the Telescope configuration parameters

	configure(configParams);

	// Set the heartbeat interval of the telescope.
	// The Telescope properties (e.g. the coordinates of the optical axis) are evolving in time, 
    // for example because of thermo-elastic drift, or because of the jitter of the platform it 
    // is mounted on. To properly track these changes one has to use a small enough timestep, 
    // which is called the "heartbeat" interval of the Telescope. Because Telescope depends on 
    // other components, like Platform which in turn may also have a certain heartbeat, the
    // 'global' heartbeat of Telescope is the minimum of its own intrinsic heartbeat and the
    // heartbeat of all the components it depends on.

	if (driftTimeScale != 0.0)
	{
		heartbeatInterval = driftTimeScale / 20.0;
	}
}










/**
 * Destructor.
 */

Telescope::~Telescope()
{

}











/**
 * \brief Configure the Telescope object using the ConfigurationParameters
 * 
 * \param configParam: the configuration parameters 
 **/

 void Telescope::configure(ConfigurationParameters &configParams)
 {
 	// Configuration parameters for the Telescope

 	lightCollectingArea     = configParams.getDouble("Telescope/LightCollectingArea");  
	transmissionEfficiency  = configParams.getDouble("Telescope/TransmissionEfficiency"); 
	FOVsolidAngle           = sqDeg2sr(configParams.getDouble("Telescope/FOVSquareDegrees"));  
	driftYawRms             = deg2rad(configParams.getDouble("Telescope/DriftYawRms") / 3600.);             
    driftPitchRms           = deg2rad(configParams.getDouble("Telescope/DriftPitchRms") / 3600.);           
    driftRollRms            = deg2rad(configParams.getDouble("Telescope/DriftRollRms") /3600.);            
    driftTimeScale          = configParams.getDouble("Telescope/DriftTimeScale");    
}










/**
 * \brief Update the telescope's pointing coordinates, by evolving the 
 *        the pointing coordinates (due to e.g. jitter or thermo-elastic variations)
 *        until time point 'time'.  
 */ 

void Telescope::updatePointingCoordinates(double time)
{
    // We can't turn back the clock, so 'time' needs to be in the future.

	if (time < internalTime)
	{
		Log.error("Telescope: cannot update pointing coordinates to time in the past");
		exit(1);
	}

    // If the give time equals exactly the current internal time, then there is nothing to update.

    if (time == internalTime)
    {
        return;
    }

    // There is currently no thermo-elastic variations in Telescope, so simply copy the 
    // pointing coordinates from platform

    tie(currentAlphaOpticalAxis, currentDeltaOpticalAxis) = platform.getPointingCoordinates(time);

    Log.info("Telescope: At time " + to_string(time) + ": (RA, dec) = (" 
                                   + to_string(rad2deg(currentAlphaOpticalAxis)) + ", " 
                                   + to_string(rad2deg(currentDeltaOpticalAxis)) + ")");

    // Update the internal clock

    internalTime = time;
	
    return;
}









/**
 * \brief Return the current values of the equatorial coordinates of the optical axis of the telescope
 * 
 * \return a pair (alphaOpticalAxis, deltaOpticalAxis)  in [rad]
 */

pair<double, double> Telescope::getCurrentPointingCoordinates()
{
	return make_pair(currentAlphaOpticalAxis, currentDeltaOpticalAxis);
}











/**
 * \brief Return the transmission efficiency of the telescope (number between 0 and 1)
 * 
 */

double Telescope::getTransmissionEfficiency()
{
	return transmissionEfficiency;
}












/**
 * \brief Return the effective light collecting area (in [m^2])
 * 
 */


double Telescope::getLightCollectingArea()
{
	return lightCollectingArea;
}










/**
 * \brief Return the solid angle covered by the FOV of 1 telescope [sr]
 */

double Telescope::getFOVsolidAngle()
{
	return FOVsolidAngle;
}











/**
 * \brief Return the equatorial sky coordinates of the optical axis of this telescope given the pointing
 *        coordinates of the (roll axis of the) platform.
 * 
 * \param alphaPlatform   Right Ascension of the pointing axis of the platform [rad]
 * \param deltaPlatform   Declination of the pointing axis of the platform     [rad]
 * 
 * \return (alphaOpticalAxis, deltaOpticalAxis)  equatorial sky coordinates of the optical axis [rad]
 */

pair<double, double> Telescope::platformToTelescopePointingCoordinates(double alphaPlatform, double deltaPlatform)
{
    // We currently assume that the telescope is perfectly aligned with the platform pointing (roll) axis
    
    const double alphaOpticalAxis = alphaPlatform;
    const double deltaOpticalAxis = deltaPlatform;

    return make_pair(alphaOpticalAxis, deltaOpticalAxis);
}













/**
 * \brief  Compute the cartesian coordinates of a point w.r.t. focal plane reference frame, given the cartesian
 *         coordinates in the spacecraft reference system.
 *         
 * \details See technical note PLATO-KUL-PL-TN-001 (De Ridder et al.)
 * 
 * \param xSC  X-coordinate in the spacecraft reference frame
 * \param ySC  Y-coordinate in the spacecraft reference frame
 * \param zSC  Z-coordinate in the spacecraft reference frame
 * \return (xFP, yFP, zFP)  Cartesian coordinates in the Focal Plane (NOT xFPprime, yFPprime, zFPprime) 
 */

tuple<double, double, double> Telescope::spacecraftToFocalPlaneCoordinates(const double xSC, const double ySC, const double zSC)
{
    // Currently the spacecraft frame equals the focal plane reference frame
    // => jitter axis = optical axis

    const double xFP = xSC;
    const double yFP = ySC;
    const double zFP = zSC;

    return make_tuple(xFP, yFP, zFP);
}











/**
 * \brief  Compute the cartesian coordinates of a point w.r.t. spacecraft reference frame, given the cartesian
 *         coordinates in the focal plane reference system.
 *         
 * \details See technical note PLATO-KUL-PL-TN-001 (De Ridder et al.)
 * 
 * \param xFP  X-coordinate in the spacecraft reference frame  (not xFPprime)
 * \param yFP  Y-coordinate in the spacecraft reference frame  (not yFPprime)
 * \param zFP  Z-coordinate in the spacecraft reference frame  (not zFPprime)
 * \return (xSC, ySC, zSC)  Cartesian coordinates in the Focal Plane reference frame.
 */

tuple<double, double, double> Telescope::focalPlaneToSpacecraftCoordinates(const double xFP, const double yFP, const double zFP)
{
    // Currently the spacecraft frame equals the focal plane reference frame
    // => jitter axis = optical axis

    const double xSC = xFP;
    const double ySC = yFP;
    const double zSC = zFP;

    return make_tuple(xSC, ySC, zSC);
}

