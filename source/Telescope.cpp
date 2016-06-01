#include "Telescope.h"

/**
 * Constructor
 * 
 * \param configurationParameters: Configuration parameters for the telescope.
 * \param hdf5File                 Output HDF5 file.
 * \param Platform:                Platform on which the telescope is mounted
 * 
 */

Telescope::Telescope(ConfigurationParameters &configParams, HDF5File &hdf5File, Platform &platform, DriftGenerator &driftGenerator)
: HDF5Writer(hdf5File), originalAzimuthAngle(0.0), originalTiltAngle(0.0), currentAzimuthAngle(0.0), currentTiltAngle(0.0),
  useDrift(true), originalFocalPlaneOrientation(0.0), currentFocalPlaneOrientation(0.0), internalTime(0.0), 
  driftGenerator(driftGenerator), platform(platform)
{
    // Initialise the HDF5 group(s) in the output file

    initHDF5Groups();

	// Retrieve the Telescope configuration parameters

	configure(configParams);

	// Initialise the heartbeat interval of the telescope.
	// The Telescope properties (e.g. the coordinates of the optical axis) are evolving in time, 
    // for example because of thermo-elastic drift, or because of the jitter of the platform it 
    // is mounted on. To properly track these changes one has to use a small enough timestep, 
    // which is called the "heartbeat" interval of the Telescope. Because Telescope depends on 
    // other components, like Platform which in turn may also have a certain heartbeat, the
    // 'global' heartbeat of Telescope is the minimum of its own intrinsic heartbeat and the
    // heartbeat of all the components it depends on.

    heartbeatInterval = min(platform.getHeartbeatInterval(), driftGenerator.getHeartbeatInterval());


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


    originalAzimuthAngle    = deg2rad(configParams.getDouble("Telescope/AzimuthAngle"));           // [rad]
    originalTiltAngle       = deg2rad(configParams.getDouble("Telescope/TiltAngle"));              // [rad]
 	lightCollectingArea     = configParams.getDouble("Telescope/LightCollectingArea") * 1.e-4;     // [m^2]  
	transmissionEfficiency  = configParams.getDouble("Telescope/TransmissionEfficiency");          // [unitless]
    useDrift                = configParams.getBoolean("Telescope/UseDrift");

    currentAzimuthAngle = originalAzimuthAngle;
    currentTiltAngle = originalTiltAngle;

    // The focal plane orientation gamma_FP used to be in the Camera class, but got moved out because 
    // it corresponds to the telescope roll orientation, which is susceptible to thermo-elastic drift, 
    // and which only Telescope can vary. Camera can access this variable throught the method
    //   Telescope::getCurrentFocalPlaneOrientation()

    originalFocalPlaneOrientation  = deg2rad(configParams.getDouble("Camera/FocalPlaneOrientation"));  // [rad]
    currentFocalPlaneOrientation = originalFocalPlaneOrientation;
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
        hdf5File.writeArray("/Telescope/", "Time",           historyTime.data(),    historyTime.size());
        hdf5File.writeArray("/Telescope/", "TelescopeRA",    historyRA.data(),      historyRA.size());         // [deg]
        hdf5File.writeArray("/Telescope/", "TelescopeDec",   historyDec.data(),     historyDec.size());        // [deg]
        hdf5File.writeArray("/Telescope/", "TelescopeYaw",   historyYaw.data(),     historyYaw.size());        // [arcsec]
        hdf5File.writeArray("/Telescope/", "TelescopePitch", historyPitch.data(),   historyPitch.size());      // [arcsec]
        hdf5File.writeArray("/Telescope/", "TelescopeRoll",  historyRoll.data(),    historyRoll.size());       // [arcsec]
        hdf5File.writeArray("/Telescope/", "Azimuth",        historyAzimuth.data(), historyAzimuth.size());    // [deg]
        hdf5File.writeArray("/Telescope/", "Tilt",           historyTilt.data(),    historyTilt.size());       // [deg]
        hdf5File.writeArray("/Telescope/", "FocalPlaneOrientation", historyFocalPlaneOrientation.data(),  historyFocalPlaneOrientation.size());   // [deg]
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
        Log.debug("Telescope: updatePointingCoordinates: coordinates up-to-date for requested time " + to_string(time));
        Log.debug("Telescope: At time " + to_string(time) + ": (azimuth, tilt, roll orient) = (" 
                                        + to_string(rad2deg(currentAzimuthAngle)) + ", " 
                                        + to_string(rad2deg(currentTiltAngle)) + ", " 
                                        + to_string(rad2deg(currentFocalPlaneOrientation)) + ") deg");

        Log.info("Telescope: At time " + to_string(time) + ": (RA, dec) = (" 
                               + to_string(rad2deg(currentAlphaOpticalAxis)) + ", " 
                               + to_string(rad2deg(currentDeltaOpticalAxis)) + ")");
       
        // If we haven't saved the optical axis coordinates for this time point yet, do so.

        if (historyTime.empty())
        {
            historyTime.push_back(time);
            historyRA.push_back(rad2deg(currentAlphaOpticalAxis));                            // [deg]
            historyDec.push_back(rad2deg(currentDeltaOpticalAxis));                           // [deg]
            historyYaw.push_back(0.0);                                                        // [arcsec]
            historyPitch.push_back(0.0);                                                      // [arcsec]
            historyRoll.push_back(0.0);                                                       // [arcsec]
            historyAzimuth.push_back(rad2deg(currentAzimuthAngle));                           // [deg]
            historyTilt.push_back(rad2deg(currentTiltAngle));                                 // [deg]
            historyFocalPlaneOrientation.push_back(rad2deg(currentFocalPlaneOrientation));    // [deg]
        }

        return;
    }

 
    // Update the azimuth, tilt, and roll orientation of the telescope which change in time due to a thermo-elastic drift

    double yaw=0.0, pitch=0.0, roll=0.0;

    if (useDrift)
    {
        const double timeInterval = time - internalTime;
        tie(yaw, pitch, roll) = driftGenerator.getNextYawPitchRoll(timeInterval);   // [rad]

        currentAzimuthAngle = originalAzimuthAngle + yaw;
        currentTiltAngle = originalTiltAngle + pitch;
        currentFocalPlaneOrientation = originalFocalPlaneOrientation + roll;

        Log.debug("Telescope: At time " + to_string(time) + ": (yaw, pitch, roll) = (" 
                                        + to_string(rad2deg(yaw)*3600.) + ", " 
                                        + to_string(rad2deg(pitch)*3600.) + ", " 
                                        + to_string(rad2deg(roll)*3600.) + ") arcsec");
    }
    else
    {
        Log.info("Telescope: Ignoring drift, telescope (yaw, pitch, roll) = (0.0, 0.0, 0.0)");
    }


    // Log the current telescope orientation of the telescope on the platform

    Log.debug("Telescope: At time " + to_string(time) + ": (azimuth, tilt, roll orient) = (" 
                                    + to_string(rad2deg(currentAzimuthAngle)) + ", " 
                                    + to_string(rad2deg(currentTiltAngle)) + ", " 
                                    + to_string(rad2deg(currentFocalPlaneOrientation)) + ") deg");




    // The telescope's optical axis does not need to be aligned with the platform's pointing axis,
    // but is usually oriented differently. Compute the equatorial sky coordinates of the telescope's
    // optical axis.

    double platformPointingRA, platformPointingDec;
    tie(platformPointingRA, platformPointingDec) = platform.getPointingCoordinates(time);
    tie(currentAlphaOpticalAxis, currentDeltaOpticalAxis) = platformToTelescopePointingCoordinates(platformPointingRA, platformPointingDec);


    // Log the current telescope pointing coordinates

    Log.debug("Telescope: At time " + to_string(time) + ": (RA, dec) = (" 
                                   + to_string(rad2deg(currentAlphaOpticalAxis)) + ", " 
                                   + to_string(rad2deg(currentDeltaOpticalAxis)) + ")");


    // Save the pointing values to write to HDF5

    historyTime.push_back(time);
    historyRA.push_back(rad2deg(currentAlphaOpticalAxis));                            // [deg]
    historyDec.push_back(rad2deg(currentDeltaOpticalAxis));                           // [deg]
    historyYaw.push_back(rad2deg(yaw) * 3600.);                                       // [arcsec]
    historyPitch.push_back(rad2deg(pitch) * 3600.);                                   // [arcsec]
    historyRoll.push_back(rad2deg(roll) * 3600.);                                     // [arcsec]
    historyAzimuth.push_back(rad2deg(currentAzimuthAngle));                           // [deg]
    historyTilt.push_back(rad2deg(currentTiltAngle));                                 // [deg]
    historyFocalPlaneOrientation.push_back(rad2deg(currentFocalPlaneOrientation));    // [deg]


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
 * \brief Return the heartbeat interval of Telescope. It's the largest interval with which 
 *        small Telescope variations (e.g. due to Platform Jitter or due to thermo-elastic
 *        drift variations) can be tracked.
 * 
 */

double Telescope::getHeartbeatInterval()
{
    return min(platform.getHeartbeatInterval(), driftGenerator.getHeartbeatInterval());
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
 * \brief  Return the current value of the focal plane orientation. The current value
 *         differs from the original value in the input file due to thermo-elastic drift.
 * 
 * \details  This variable used to be in the Camera class, but got moved out because it 
 *           corresponds to the telescope roll orientation, which is susceptible to thermo-
 *           elastic drift, and which only Telescope can vary. 
 *            
 * \note  See also PLATO-KUL-PL-TN-0001
 * 
 * \return current value of the focal plane orientation gamma_FP  [rad]
 */

double Telescope::getCurrentFocalPlaneOrientation()
{
    return currentFocalPlaneOrientation;
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

    double cosAngle = cos(currentAzimuthAngle);
    double sinAngle = sin(currentAzimuthAngle);

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

    cosAngle = cos(currentTiltAngle);
    sinAngle = sin(currentTiltAngle);

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








