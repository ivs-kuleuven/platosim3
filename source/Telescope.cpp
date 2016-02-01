#include "Telescope.h"

/**
 * Constructor
 * 
 * \param configurationParameters: Configuration parameters for the telescope.
 * \param Platform:                Platform on which the telescope is mounted
 * \param hdf5File                 Output HDF5 file.
 * 
 */

Telescope::Telescope(ConfigurationParameters &configParams, HDF5File &hdf5File)
: HDF5Writer(hdf5File)
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

	// Also destroy the platform

}





/**
 * Method that updates the pointing coordinates of the telescope (i.e. the
 * equatorial coordinates of the optical axis) after the application of jitter
 * (of the platform on which the telescope is mounted) and the thermo-elastic
 * variations (of the telescope itself).
 *
 * @param &alphaOpticalAxis: Right ascension of the optical axis, before applying
 *        the next jitter step and taking thermo-elastic variations into account
 *        [degees?/radians?]
 *
 * @param &deltaOpticalAxis: Declination of the optical axis, before applying the
 *        next jitter step and taking thermo-elastic variations into account
 *        [degrees?/radians?]
 */
double Telescope::updatePointingCoordinates(double &alphaOpticalAxis,
		double &deltaOpticalAxis, double currentTime)
{

	// Jitter -> platform
	// Ask the platform to update the pointing coordinates (i.e. the equatorial
	// coordinates of the telescope axis) and the current time, taking the
	// jitter into account

//	this->getPlatform().updatePointingCoordinates(alphaOpticalAxis, deltaOpticalAxis, currentTime);

	// Thermo-elastic variations
	// Update the coordinates of the displaced (because of the jitter) optical
	// axis, taking the thermo-elastic variations into account

//	this->getPlatform().updatePointingCoordinates(alphaOpticalAxis, deltaOpticalAxis, currentTime);

	return currentTime;

}





// Platform

/**
 * Method that mounts the telescope on the given platform.
 *
 * @param platform: Platform on which to mount the telescope.
 * @type platform: Platform
 */
void Telescope::setPlatform(Platform platform) {

	this->platform = platform;

}





/**
 * Method that returns the platform on which the telescope is mounted.
 *
 * @return Platform on which the telescope is mounted.
 * @rtype Platform
 */
Platform Telescope::getPlatform()
{

	return this->platform;

}

// Light collecting area

/**
 * Method that sets the light collecting area of the telescope to the given area.
 *
 * @param lightCollectingArea: Light collecting area to use for the telescope [cm^2].
 * @type lightCollectingArea: double
 */
void Telescope::setLightCollectingArea(double lightCollectingArea)
{

	this->lightCollectingArea = lightCollectingArea;

}





/**
 * Method that returns the light collecting area of the telescope [cm^2].
 *
 * @return Light collecting area of the telescope [cm^2].
 * @rtype double
 */
unsigned double Telescope::getLightCollectingArea()
{

	return this->lightCollectingArea;
}
