
#include "Platform.h"



/**
 * \brief Constructor
 * 
 * \param configParams        Configuration parameters as read from the (e.g. yaml) inputfile
 * \param hdf5File            HDF5 file where to write any output
 * \param jitterGenerator     Generator of the yaw, pitch roll variations due to spacecraft jitter
 * 
 * \note  The jitterGenerator has been configured in Simulation::Simulation()
 */

Platform::Platform(ConfigurationParameters configParams, HDF5File &hdf5File, JitterGenerator &jitterGenerator)
: HDF5Writer(hdf5File), useJitter(true), internalTime(0.0), jitterGenerator(jitterGenerator)
{
    // Initialise the HDF5 group(s) in the output file

    initHDF5Groups();

    // Configure the Platfrom object

    configure(configParams);
}







/**
 * \brief Destructor
 */

Platform::~Platform()
{
    flushOutput();
}








/**
 * \brief Configure the Platform object
 * 
 * \param configParams  The configuration parameters from the input parameters file
 */

void Platform::configure(ConfigurationParameters &configParams)
{
    useJitter   = configParams.getBoolean("Platform/UseJitter");
    originalRA  = deg2rad(configParams.getDouble("ObservingParameters/RApointing"));            
    originalDec = deg2rad(configParams.getDouble("ObservingParameters/DecPointing"));
    raSun       = deg2rad(configParams.getDouble("ObservingParameters/RASun"));            
    decSun      = deg2rad(configParams.getDouble("ObservingParameters/DecSun"));
         
    currentRA   = originalRA;
    currentDec  = originalDec;
}












/**
 * \brief Creates the group(s) in the HDF5 file where the ACS information will be stored. 
 *        These group(s) have to be created once, at the very beginning.
 */

void Platform::initHDF5Groups()
{
    Log.debug("Platform: initialising HDF5 groups");

    hdf5File.createGroup("/ACS");
}











/**
 * \brief Write all recorded information to the HDF5 output file
 */

void Platform::flushOutput()
{
    Log.info("Platform: Flushing output to HDf5 file.");

    if ( ! hdf5File.hasGroup("ACS") )
    {
        Log.warning("Platform.flushOutput: HDF5 file has no ACS group, cannot flush Platform information.");
        return;
    }
    

     if (!historyTime.empty())
     {
        hdf5File.writeArray("/ACS/", "Time",        historyTime.data(),  historyTime.size());
        hdf5File.writeArray("/ACS/", "PlatformRA",  historyRA.data(),    historyRA.size());
        hdf5File.writeArray("/ACS/", "PlatformDec", historyDec.data(),   historyDec.size());
        hdf5File.writeArray("/ACS/", "Yaw",         historyYaw.data(),   historyYaw.size());
        hdf5File.writeArray("/ACS/", "Pitch",       historyPitch.data(), historyPitch.size());
        hdf5File.writeArray("/ACS/", "Roll",        historyRoll.data(),  historyRoll.size());
     }
     else
     {
        Log.warning("Platform: No ACS history to flush to HDF5 file.");
     }
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
    originalRA = rightAscencsion / unit;
    originalDec = declination / unit;

    currentRA  = originalRA;
    currentDec = originalDec;
}












/**
 * \brief Update the pointing coordinates (of the roll axis) of the spacecraft
 * 
 * \param time [s]
 *  
 * \return Nothing
 */

void Platform::updatePointingCoordinates(double time)
{
    double yaw=0.0, pitch=0.0, roll=0.0;

    if (useJitter)
    {
        // Check if the request is going backwards in time. If so: complain.

        double timeInterval = time - internalTime;

        if (timeInterval < 0.0)
        {
            Log.warning("Platform: getPointingCoordinates() at time before previous request: Not Implemented. ");
            return;
        }


        // Check if the given 'time' is the same as the last one we processed (and kept).
        // If so, nothing has to change.

        if (!historyTime.empty())
        {
            if (time == historyTime.back())
            {
                Log.debug("Platform: getPointingCoordinates: coordinates up-to-date for requested time " + to_string(time));
                Log.debug("Platform: At time " + to_string(time) + ": (RA, dec) = (" 
                                               + to_string(historyRA.back()) + ", " 
                                               + to_string(historyDec.back()) + ")");

                return;
            }
        }

        // We're know in the case that we haven't processed the given time point yet.
        // Let the platfrom jitter until 'time'. Yaw, pitch, and roll are in [rad]

        tie(yaw, pitch, roll) = jitterGenerator.getNextYawPitchRoll(timeInterval);

        Log.debug("Platform: At time " + to_string(time) + ": (yaw, pitch, roll) = (" 
                                       + to_string(rad2deg(yaw)*3600.) + ", " 
                                       + to_string(rad2deg(pitch)*3600.) + ", " 
                                       + to_string(rad2deg(roll)*3600.) + ") arcsec");


        // The roll axis (= unit vector in z-direction in SC reference frame) will have slightly 
        // rotated due to jitter. Find out the cartesian coordinates of the _new_ jitter axis in 
        // the SpaceCraft reference frame of the original pointing. 

        arma::colvec zUnitBeforeJitter = {0.0, 0.0, 1.0};
        arma::colvec zUnitAfterJitter = rotateYawPitchRoll(zUnitBeforeJitter, yaw, pitch, roll);

        // Compute the celestial equatorial cartesian coordinates of the new roll axis
        // This requires the original pointing coordinates of the platform.

        const bool useOriginalPointingCoordinates = true;
        const arma::colvec zUnitAfterJitterEQ = spacecraftToEquatorialCoordinates(zUnitAfterJitter, useOriginalPointingCoordinates);

        // Convert from cartesian to celestial equatorial coordinates

        const double x = zUnitAfterJitterEQ(0);
        const double y = zUnitAfterJitterEQ(1);
        const double z = zUnitAfterJitterEQ(2);

        // Only now update the internal platform pointing coordinates to the ones after the jitter step
        // Note: r should 1.0, as rotations don't change the length of the unit vector

        const double r = sqrt(x*x+y*y+z*z);
        currentDec = PI / 2.0 - acos(z/r);
        currentRA = atan2(y, x);
        if (currentRA < 0.0) currentRA += 2 * PI; 
    }
    else
    {
        Log.info("Platform: Ignoring jitter, platform (yaw, pitch, roll) = (0.0, 0.0, 0.0)");
        yaw = 0.0;
        pitch = 0.0;
        roll = 0.0;
        currentRA = originalRA;
        currentDec = originalDec;
    }

    Log.debug("Platform: At time " + to_string(time) + ": (RA, dec) = (" 
                                   + to_string(rad2deg(currentRA)) + ", " 
                                   + to_string(rad2deg(currentDec)) + ")");

    // Update the internal clock

    internalTime = time;

    // Store the computed values so that they can later be saved to HDF5
    // RA & Dec are saved in degrees. Yaw, pitch, roll in arcsec.

    historyTime.push_back(time);
    historyRA.push_back(rad2deg(currentRA));
    historyDec.push_back(rad2deg(currentDec));
    historyYaw.push_back(rad2deg(yaw) * 3600.);
    historyPitch.push_back(rad2deg(pitch) * 3600.);
    historyRoll.push_back(rad2deg(roll) * 3600.);

    // That's it

    return;
}













/**
 * \brief Return the current (RA, Dec) pointing coordinates of the spacecraft
 * 
 * \return (RA, Dec)  Right Ascension and declination [rad]
 */

pair<double, double> Platform::getCurrentPointingCoordinates()
{
    return make_pair(currentRA, currentDec);
}













/**
 * \brief Return the current (RA, Dec) pointing coordinates of the spacecraft
 * 
 * \return (RA, Dec)  Right Ascension and declination [rad]
 */

pair<double, double> Platform::getInitialPointingCoordinates()
{
    return make_pair(originalRA, originalDec);
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














/**
 * \brief Compute 3D cartesian coordinates in the celestial equatorial reference frame,
 *        given the 3D cartesian coordinates in the spacecraft (SC) reference frame
 * 
 * \param vecSC   (xSC, ySC, zSC): cartesian coordinates of the point in the spacecraft reference frame
 * \param useOriginalPointingCoordinates  If true: use original pointing coordinates (before jitter started)
 *                                        If false: use current pointing coordinates (affected by jitter)
 *
 * \return coordEQ  (xEQ, yEQ, zEQ): cartesian coordinates in the celestial equatorial reference frame
 */

arma::colvec Platform::spacecraftToEquatorialCoordinates(arma::colvec &vecSC, bool useOriginalPointingCoordinates)
{

    // Follow the user on which pointing coordinates to use: the current or the original ones

    double RA, dec;

    if (useOriginalPointingCoordinates)
    {
        RA = originalRA;
        dec = originalDec;
    }
    else
    {
        RA = currentRA;
        dec = currentDec;
    }

    // Compute the equatorial coordinates of each of the unit vectors corresponding to the X, Y, and Z axis
    // of the spacecraft reference frame. The z-axis is pointing towards the targets, the x-axis points towards
    // the highest point of the sun shield which is pointing towards the Sun.

    double deltax = atan(- cos(RA-raSun) / tan(dec));

    arma::colvec unitzSC = {cos(dec)*cos(RA), cos(dec)*sin(RA), sin(dec)};
    arma::colvec unitxSC = {cos(deltax)*cos(raSun), cos(deltax)*sin(raSun), sin(deltax)};
    arma::colvec unitySC = arma::cross(unitzSC, unitxSC);

    // Compute the rotation matrix to convert cartesian coordinates in the equatorial reference frame to 
    // cartesian coordinates in the spacecraft reference frame

    arma::mat rotSC2EQ;
    rotSC2EQ << unitxSC[0] << unitySC[0] << unitzSC[0] << arma::endr
             << unitxSC[1] << unitySC[1] << unitzSC[1] << arma::endr
             << unitxSC[2] << unitySC[2] << unitzSC[2] << arma::endr;
    

    // Transform the cartesian spacecraft coordinates to cartesian equatorial coordinates

    arma::colvec vecEQ = rotSC2EQ * vecSC;

    // That's it

    return vecEQ;
}












/**
 * \brief Transforms the coordinates of a vector after yaw, pitch, and roll rotations.
 * 
 * \details We define the yaw, the pitch, and the roll as rotation angles around respectively 
 *          the X_SC, Y_SC, and Z_SC axes, (SC standing for spacecraft reference frame) such that 
 *          the angles increase with a clockwise rotation, when looking along the positive axes. 
 *          First a roll rotation is done around the Z_SC axis, then a pitch rotation is done 
 *          around the rotated Y_SC axis, and finally a yaw rotation is done around the already
 *          twice-rotated X_SC axis. The combined rotation matrix is: 
 *                \f[ R(yaw, pitch, roll) = R(yaw) * R(pitch) * R(roll) \f]
 *          Which transforms
 *                \f[ vector\_after = R(yaw, pitch, roll) * vector\_before \f]
 *          where vector_before and vector_after are the coordinates of the same point before and
 *          after the yaw.pitch,roll rotations. 
 *          
 * \note See also Technical Note PLATO-KUL-PL-TN-001
 * 
 * \param coord  Coordinates of the vector in the spacecraft reference frame before the rotation 
 * \param yaw    Yaw angle   [rad]
 * \param pitch  Pitch angle [rad]
 * \param roll   Roll angle  [rad]
 * \return       Coordinates of the vector in the (old) spacecraft reference frame after the rotations
 */

arma::colvec Platform::rotateYawPitchRoll(arma::colvec coord, const double yaw, const double pitch, const double roll)
{
    // Some handy abbreviations

    const double sinYaw   = sin(yaw);
    const double cosYaw   = cos(yaw);
    const double sinPitch = sin(pitch);
    const double cosPitch = cos(pitch);
    const double sinRoll  = sin(roll);
    const double cosRoll  = cos(roll);

    // The rotation matrices

    arma::mat Ryaw;
    Ryaw << 1.0  <<  0.0    <<     0.0   << arma::endr
         << 0.0  <<  cosYaw <<  -sinYaw  << arma::endr
         << 0.0  <<  sinYaw <<   cosYaw  << arma::endr;

    arma::mat Rpitch;
    Rpitch <<  cosPitch  <<  0.0  <<  sinPitch  << arma::endr
           <<   0.0      <<  1.0  <<    0.0     << arma::endr
           << -sinPitch  <<  0.0  <<  cosPitch  << arma::endr;

    arma::mat Rroll;
    Rroll << cosRoll  <<  -sinRoll  <<  0.0  << arma::endr
          << sinRoll  <<   cosRoll  <<  0.0  << arma::endr
          <<  0.0     <<     0.0    <<  1.0  << arma::endr; 

    // Do the transformation

    return Ryaw  * Rpitch * Rroll * coord;
}

