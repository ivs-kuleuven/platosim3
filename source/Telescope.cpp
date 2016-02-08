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

}











/**
 * \brief Configure the Telescope object using the ConfigurationParameters
 * 
 * \param configParam: the configuration parameters 
 **/

 void Telescope::configure(ConfigurationParameters &configParam)
 {
 	// Configuration parameters for the Telescope

 	alphaOpticalAxis        = deg2rad(configParam.getDouble("ObservingParameters/RApointing"));
	deltaOpticalAxis        = deg2rad(configParam.getDouble("ObservingParameters/DecPointing"));        
	lightCollectingArea     = configParam.getDouble("Telescope/lightCollectingArea");     
	transmissionEfficiency  = configParam.getDouble("Telescope/TransmissionEfficiency");  
	driftYawRms             = configParam.getDouble("Telescope/DriftYawRms");             
    driftPitchRms           = configParam.getDouble("Telescope/DriftPitchRms");           
    driftRollRms            = configParam.getDouble("Telescope/DriftRollRms");            
    driftTimeScale          = configParam.getDouble("Telescope/DriftTimeScale");    
}










/**
 * \brief Update the telescope's pointing coordinates. Over the given 'timeInterval' they may
 *        change due to platform jitter or thermo-elastic variations.
 *          
 */ 

void Telescope::updatePointingCoordinates(double timeInterval)
{
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

