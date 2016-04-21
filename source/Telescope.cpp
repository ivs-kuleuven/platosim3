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
: HDF5Writer(hdf5File), azimuthAngle(0.0), tiltAngle(0.0), internalTime(0.0), platform(platform)
{
    // Initialise the HDF5 group(s) in the output file

    initHDF5Groups();

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

    // Initialize the current position of the optical axis of the telescope given the 
    // platform pointing axis.

    double alphaPlatform, deltaPlatform;
    tie(alphaPlatform, deltaPlatform) = platform.getPointingCoordinates(internalTime);
    tie(currentAlphaOpticalAxis, currentDeltaOpticalAxis) = platformToTelescopePointingCoordinates(alphaPlatform, deltaPlatform);
}










/**
 * Destructor.
 */

Telescope::~Telescope()
{
    flushOutput();
}











/**
 * \brief Configure the Telescope object using the ConfigurationParameters
 * 
 * \param configParam: the configuration parameters 
 **/

 void Telescope::configure(ConfigurationParameters &configParams)
 {
 	// Configuration parameters for the Telescope


    azimuthAngle            = deg2rad(configParams.getDouble("Telescope/AzimuthAngle"));           // [rad]
    tiltAngle               = deg2rad(configParams.getDouble("Telescope/TiltAngle"));              // [rad]
 	lightCollectingArea     = configParams.getDouble("Telescope/LightCollectingArea") * 1.e-4;     // [m^2]  
	transmissionEfficiency  = configParams.getDouble("Telescope/TransmissionEfficiency");          // [unitless]
	FOVsolidAngle           = sqDeg2sr(configParams.getDouble("Telescope/FOVSquareDegrees"));      // [sr]
	driftYawRms             = deg2rad(configParams.getDouble("Telescope/DriftYawRms") / 3600.);    // [rad]         
    driftPitchRms           = deg2rad(configParams.getDouble("Telescope/DriftPitchRms") / 3600.);  // [rad]         
    driftRollRms            = deg2rad(configParams.getDouble("Telescope/DriftRollRms") /3600.);    // [s]        
    driftTimeScale          = configParams.getDouble("Telescope/DriftTimeScale");    
}








/**
 * \brief Creates the group(s) in the HDF5 file where the Telescope pointing information
 *        will be stored. These group(s) have to be created once, at the very beginning.
 */

void Telescope::initHDF5Groups()
{
    Log.debug("Telescope: initialising HDF5 groups");

    hdf5File.createGroup("/Telescope");
}









/**
 * \brief Write all recorded information to the HDF5 output file
 */

void Telescope::flushOutput()
{
    Log.info("Telescope: Flushing output to HDf5 file.");

    if ( ! hdf5File.hasGroup("Telescope") )
    {
        Log.warning("Telescope::flushOutput: HDF5 file has no Telescope group, cannot flush Telescope information.");
        return;
    }
    

     if (!historyTime.empty())
     {
        hdf5File.writeArray("/Telescope/", "Time",         historyTime.data(),  historyTime.size());
        hdf5File.writeArray("/Telescope/", "TelescopeRA",  historyRA.data(),    historyRA.size());
        hdf5File.writeArray("/Telescope/", "TelescopeDec", historyDec.data(),   historyDec.size());
     }
     else
     {
        Log.warning("Telescope: No telescope pointing history to flush to HDF5 file.");
     }
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
        Log.info("Telescope: At time " + to_string(time) + ": (RA, dec) = (" 
                               + to_string(rad2deg(currentAlphaOpticalAxis)) + ", " 
                               + to_string(rad2deg(currentDeltaOpticalAxis)) + ")");
       
        // If we haven't saved the optical axis coordinates for this time point yet, do so.

        if (historyTime.empty())
        {
            historyTime.push_back(time);
            historyRA.push_back(rad2deg(currentAlphaOpticalAxis));
            historyDec.push_back(rad2deg(currentDeltaOpticalAxis));
        }
        else
        {   
            if (historyTime.back() != time)
            {
                historyTime.push_back(time);
                historyRA.push_back(rad2deg(currentAlphaOpticalAxis));
                historyDec.push_back(rad2deg(currentDeltaOpticalAxis));
            }
        }

        return;
    }

    // Get the updated pointing coordinates of the platform

    double platformPointingRA, platformPointingDec;
    tie(platformPointingRA, platformPointingDec) = platform.getPointingCoordinates(time);

    // The telescope's optical axis does not need to be aligned with the platform's pointing axis,
    // but is usually oriented differently. Compute the equatorial sky coordinates of the telescope's
    // optical axis.

    tie(currentAlphaOpticalAxis, currentDeltaOpticalAxis) = platformToTelescopePointingCoordinates(platformPointingRA, platformPointingDec);

    // Apply a thermo-elastic drift 

    // TBD

    // Log the current telescope pointing coordinates

    Log.info("Telescope: At time " + to_string(time) + ": (RA, dec) = (" 
                                   + to_string(rad2deg(currentAlphaOpticalAxis)) + ", " 
                                   + to_string(rad2deg(currentDeltaOpticalAxis)) + ")");


    // Save the pointing values to write to HDF5

    historyTime.push_back(time);
    historyRA.push_back(rad2deg(currentAlphaOpticalAxis));
    historyDec.push_back(rad2deg(currentDeltaOpticalAxis));

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
 * \note  See also PLATO-KUL-PL-TN-001 for the definitions of the reference frames.
 * 
 * \param alphaPlatform   Right Ascension of the pointing axis of the platform [rad]
 * \param deltaPlatform   Declination of the pointing axis of the platform     [rad]
 * 
 * \return (alphaOpticalAxis, deltaOpticalAxis)  equatorial sky coordinates of the optical axis [rad]
 */

pair<double, double> Telescope::platformToTelescopePointingCoordinates(double alphaPlatform, double deltaPlatform)
{
    // Specify the coordinates in the equatorial reference frame of the unit vector zSC,
    // corresponding to the jitter (Z) axis of the SpaceCraft

    arma::vec zSC = {cos(alphaPlatform)*cos(deltaPlatform), sin(alphaPlatform)*cos(deltaPlatform), sin(deltaPlatform)};

    // Construct the rotation matrix for a rotation around the zSC axis over the azimuth angle
    // {ux, uy, uz} is short-hand notation

    double ux = zSC(0);
    double uy = zSC(1);
    double uz = zSC(2);

    double cosAngle = cos(azimuthAngle);
    double sinAngle = sin(azimuthAngle);

    arma::mat rotAzimuth = {{cosAngle+ux*ux*(1-cosAngle),    ux*uy*(1-cosAngle)-uz*sinAngle, ux*uz*(1-cosAngle)+uy*sinAngle},
                            {uy*ux*(1-cosAngle)+uz*sinAngle, cosAngle+uy*uy*(1-cosAngle),    uy*uz*(1-cosAngle)-ux*sinAngle},
                            {uz*ux*(1-cosAngle)-uy*sinAngle, uz*uy*(1-cosAngle)+ux*sinAngle, cosAngle+uz*uz*(1-cosAngle)}};


    // The goal of the rotZ rotation matrix is to rotate the ySC unit vector (corresponding to the y-axis in 
    // the spacecraft reference frame) in the azimuth direction of the telescope. Rather than using ySC, we use
    // another reference vector yRef as defined below. y0 is perpendicular to zSC, and has the advantage that the 
    // exact orientation of the spacecraft (i.e. in which direction the sunshield is pointing) is not needed.

    arma::vec yRef = {-sin(alphaPlatform), cos(alphaPlatform), 0.0};

    // Rotate this reference vector

    arma::vec yAzimuth = rotAzimuth * yRef;

    // Next, construct the rotation matrix for a rotation around the yAzimuth vector over the tilt angle
    // of the telescope. The tilt angle is the angle between the optical axis of the telescope and the
    // the jitter Z-axis of the platform.

    ux = yAzimuth(0);
    uy = yAzimuth(1);
    uz = yAzimuth(2);

    cosAngle = cos(tiltAngle);
    sinAngle = sin(tiltAngle);

    arma::mat rotTilt = {{cosAngle+ux*ux*(1-cosAngle),    ux*uy*(1-cosAngle)-uz*sinAngle, ux*uz*(1-cosAngle)+uy*sinAngle},
                         {uy*ux*(1-cosAngle)+uz*sinAngle, cosAngle+uy*uy*(1-cosAngle),    uy*uz*(1-cosAngle)-ux*sinAngle},
                         {uz*ux*(1-cosAngle)-uy*sinAngle, uz*uy*(1-cosAngle)+ux*sinAngle, cosAngle+uz*uz*(1-cosAngle)}};


    // Compute the unit vector zOA in the direction of the telescope's optical axis

    arma::vec zOA = rotTilt * zSC;

    // zOA now contains the cartesian coordinates of the optical axis in the equatorial reference frame. 
    // Compute the equatorial sky coordinates [rad] from the cartesian coordinates.

    const double norm = sqrt(zOA(0)*zOA(0) + zOA(1)*zOA(1) + zOA(2)*zOA(2));    // should be 1.0
    double deltaOpticalAxis = Constants::PI/2.0 - acos(zOA(2)/norm);
    double alphaOpticalAxis = atan2(zOA(1), zOA(0));

    if (alphaOpticalAxis < 0.0)
    {
        alphaOpticalAxis += 2.0 * Constants::PI;
    }

    // That's it!

    return make_pair(alphaOpticalAxis, deltaOpticalAxis);
}








