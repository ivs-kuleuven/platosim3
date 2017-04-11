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
    tie(alphaPlatform, deltaPlatform) = platform.getInitialPointingCoordinates();
    tie(currentAlphaOpticalAxis, currentDeltaOpticalAxis) = platformToTelescopePointingCoordinates(alphaPlatform, deltaPlatform);

    // Get the equatorial sky coordinates of the Sun, which is know by platform since it's pointing its sunshield towards it.

    tie(raSun, decSun) = platform.getRADecSun();
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


    originalAzimuthAngle      = deg2rad(configParams.getDouble("Telescope/AzimuthAngle"));                  // [rad]
    originalTiltAngle         = deg2rad(configParams.getDouble("Telescope/TiltAngle"));                     // [rad]
    lightCollectingArea       = configParams.getDouble("Telescope/LightCollectingArea") * 1.e-4;            // [m^2]  
    transmissionEfficiencyBOL = configParams.getDouble("Telescope/TransmissionEfficiency/BOL");
    transmissionEfficiencyEOL = configParams.getDouble("Telescope/TransmissionEfficiency/EOL");
    missionDuration           = configParams.getDouble("ObservingParameters/MissionDuration") * 31536000.0; // [s]
    useDrift                  = configParams.getBoolean("Telescope/UseDrift");

    currentAzimuthAngle = originalAzimuthAngle;
    currentTiltAngle    = originalTiltAngle;

    // The focal plane orientation gamma_FP used to be in the Camera class, but got moved out because 
    // it corresponds to the telescope roll orientation, which is susceptible to thermo-elastic drift, 
    // and which only Telescope can vary. Camera can access this variable throught the method
    //   Telescope::getCurrentFocalPlaneOrientation()

    originalFocalPlaneOrientation = deg2rad(configParams.getDouble("Camera/FocalPlaneOrientation"));  // [rad]
    currentFocalPlaneOrientation  = originalFocalPlaneOrientation;
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
 * \brief Update the telescope's pointing coordinates, by evolving the pointing coordinates 
 *        (due to e.g. jitter or thermo-elastic variations) until time point 'time'.  
 */ 



void Telescope::updateTelescopeOrientation(double time)
{
    // We can't turn back the clock, so 'time' needs to be in the future.

    if (time < internalTime)
    {
        Log.error("Telescope: cannot update telescope orientation to time in the past");
        exit(1);
    }

    // Check if the given 'time' is the same as the last one we processed (and kept).
    // If so, nothing has to change.

    if (!historyTime.empty())
    {
        if (time == historyTime.back())
        {
            Log.debug("Telescope: updateTelescopeOrientation: coordinates up-to-date for requested time " + to_string(time));
            Log.debug("Telescope: At time " + to_string(time) + ": (azimuth, tilt, roll orient) = (" 
                                            + to_string(rad2deg(currentAzimuthAngle)) + ", " 
                                            + to_string(rad2deg(currentTiltAngle)) + ", " 
                                            + to_string(rad2deg(currentFocalPlaneOrientation)) + ") deg");

            Log.info("Telescope: At time " + to_string(time) + ": (RA, dec) = (" 
                                   + to_string(rad2deg(currentAlphaOpticalAxis)) + ", " 
                                   + to_string(rad2deg(currentDeltaOpticalAxis)) + ")");

            return;
        }
    }

 
    // Update the azimuth, tilt, and roll orientation of the telescope which change in time due to a thermo-elastic drift

    double yaw=0.0, pitch=0.0, roll=0.0;

    if (useDrift)
    {
        tie(yaw, pitch, roll) = driftGenerator.getNextYawPitchRoll(time);   // [rad]

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

    double raPlatform, decPlatform;
    platform.updatePointingCoordinates(time);
    tie(raPlatform, decPlatform) = platform.getCurrentPointingCoordinates();
    tie(currentAlphaOpticalAxis, currentDeltaOpticalAxis) = platformToTelescopePointingCoordinates(raPlatform, decPlatform);


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
 * \brief Return the current azimuth and title angles of the telescope w.r.t. the Platform,
 *        the focal plane orientation angle, and the current platform pointing coordinates.
 *        Also the (RA, Dec) of the platform is given, so that all information is provided
 *        to (re)construct the telescope reference frame.
 *        
 * \details Historically the focal plane orientation was a feature of the Camera class, but since
 *          it is also serves as the roll angle of the telescope, it's included here as well.
 *          
 * \return azimuth angle of telescope [rad]
 *         tilt angle of telescope    [rad]
 *         focal plane angle          [rad]
 *         RA of platform             [rad]
 *         declination of platform    [rad]
 */

tuple<double, double, double, double, double> Telescope::getCurrentTelescopeOrientation()
{
    double raPlatform, decPlatform;
    tie(raPlatform, decPlatform) = platform.getCurrentPointingCoordinates();
    return make_tuple(currentAzimuthAngle, currentTiltAngle, currentFocalPlaneOrientation, raPlatform, decPlatform);
}

        











/**
 * \brief Return the initial azimuth and title angles of the telescope w.r.t. the Platform
 *        and the initial platform pointing coordinates, as they were specified in the configuration file.
 *        
 * \details Historically the focal plane orientation was a feature of the Camera class, but since
 *          it is also serves as the roll angle of the telescope, it's included here as well.       
 *        
 * \return azimuth angle of telescope [rad]
 *         tilt angle of telescope    [rad]
 *         focal plane angle          [rad]       
 *         RA of platform             [rad]
 *         declination of platform    [rad]
 */

tuple<double, double, double, double, double> Telescope::getInitialTelescopeOrientation()
{
    double originalRAPlatform, originalDecPlatform;
    tie(originalRAPlatform, originalDecPlatform) = platform.getInitialPointingCoordinates();
    return make_tuple(originalAzimuthAngle, originalTiltAngle, originalFocalPlaneOrientation, originalRAPlatform, originalDecPlatform);
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

double Telescope::getTransmissionEfficiency(double time)
{
    return transmissionEfficiencyBOL - (transmissionEfficiencyBOL - transmissionEfficiencyEOL) / missionDuration * time;
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
 * \brief Return the equatorial sky coordinates of the optical axis of this telescope given the pointing
 *        coordinates of the (roll axis of the) platform.
 *        
 * \note  See also PLATO-KUL-PL-TN-001 for the definitions of the reference frames.
 * 
 * \param raPlatform   Right Ascension of the pointing axis of the platform [rad]
 * \param decPlatform   Declination of the pointing axis of the platform     [rad]
 * 
 * \return (alphaOpticalAxis, deltaOpticalAxis)  equatorial sky coordinates of the optical axis [rad]
 */

pair<double, double> Telescope::platformToTelescopePointingCoordinates(double raPlatform, double decPlatform)
{
    // Compute the rotation matrix to convert cartesian coordinates in the telescope reference frame to 
    // cartesian coordinates in the spacecraft reference frame

    arma::mat rotAzimuth;
    rotAzimuth <<  cos(currentAzimuthAngle) << sin(currentAzimuthAngle) << 0 << arma::endr
               << -sin(currentAzimuthAngle) << cos(currentAzimuthAngle) << 0 << arma::endr
               <<           0               <<            0             << 1 << arma::endr;


    arma::mat rotTilt;
    rotTilt <<  cos(currentTiltAngle) << 0 << sin(currentTiltAngle) << arma::endr
            <<          0             << 1 <<            0          << arma::endr
            << -sin(currentTiltAngle) << 0 << cos(currentTiltAngle) << arma::endr;

    arma::mat rotTL2SC = rotAzimuth * rotTilt;

    // Compute the equatorial cartesian coordinates of the unit vector along the z-axis (= roll = pointing axis) of the platform.
    // The x-axis of the platform points to the highest point fof the sunshield, which is pointing to the (average) sky position
    // of the Sun.

    double deltax = atan(- cos(raPlatform-raSun) / tan(decPlatform));
    arma::colvec zSC = {cos(decPlatform)*cos(raPlatform), cos(decPlatform)*sin(raPlatform), sin(decPlatform)};
    arma::colvec xSC = {cos(deltax)*cos(raSun), cos(deltax)*sin(raSun), sin(deltax)};
    arma::colvec ySC = arma::cross(zSC, xSC);

    // Compute the rotation matrix to convert cartesian coordinates in the equatorial reference frame to 
    // cartesian coordinates in the spacecraft reference frame

    arma::mat rotSC2EQ;
    rotSC2EQ << xSC[0] << ySC[0] << zSC[0] << arma::endr
             << xSC[1] << ySC[1] << zSC[1] << arma::endr
             << xSC[2] << ySC[2] << zSC[2] << arma::endr;
    
    // Combine all the rotation matrices

    arma::mat rotTL2EQ = rotSC2EQ * rotTL2SC;

    // In the telescope reference frame, the optical axis is simply the z-axis = (0,0,1)

    arma::colvec vecTL = {0.0, 0.0, 1.0};

    // Get the equatorial coordinates of this optical axis vector.

    arma::colvec vecEQ = rotTL2EQ * vecTL;

    // Convert the cartesian equatorial coordinates to equatorial sky coordinates

    double norm = sqrt(vecEQ[0]*vecEQ[0] + vecEQ[1]*vecEQ[1] + vecEQ[2]*vecEQ[2]);
    double decOpticalAxis = Constants::PI/2.0 - acos(vecEQ[2]/norm);
    double raOpticalAxis = atan2(vecEQ[1], vecEQ[0]);

    // Ensure that the right ascension is positive

    if (raOpticalAxis < 0.0)
    {
        raOpticalAxis += 2.* Constants::PI;
    }

    // That's it!

    return make_pair(raOpticalAxis, decOpticalAxis);
}

