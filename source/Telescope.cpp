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
	driftYawRms             = configParams.getDouble("Telescope/DriftYawRms");             
    driftPitchRms           = configParams.getDouble("Telescope/DriftPitchRms");           
    driftRollRms            = configParams.getDouble("Telescope/DriftRollRms");            
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

    // Telescope depends on Platform (and its jitter) to get new pointing coordinates.
    // So first update platform.

    platform.updatePointingCoordinates(time);

    // There is currently no thermo-elastic variations in Telescope, so simply copy the 
    // pointing coordinates from platform

    tie(alphaOpticalAxis, deltaOpticalAxis) = platform.getPointingCoordinates();

    // Update the internal clock

    internalTime = time;
	
    return;
}









/**
 * \brief Return the current values of the equatorial coordinates of the optical axis of the telescope
 * 
 * \return a pair (alphaOpticalAxis, deltaOpticalAxis)  in [rad]
 */

pair<double, double> Telescope::getPointingCoordinates()
{
	return make_pair(alphaOpticalAxis, deltaOpticalAxis);
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


