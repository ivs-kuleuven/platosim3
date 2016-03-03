
#include "Platform.h"



/**
 * \brief Constructor
 * 
 * \param configParams        Configuration parameters as read from the (e.g. yaml) inputfile
 * \param hdf5File            HDF5 file where to write any output
 * \param jitterGenerator     Generator of the yaw, pitch roll variations due to spacecraft jitter
 */

Platform::Platform(ConfigurationParameters configParams, HDF5File &hdf5File, JitterGenerator &jitterGenerator)
: HDF5Writer(hdf5File), internalTime(0.0), jitterGenerator(jitterGenerator)
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
    originalRA  = deg2rad(configParams.getDouble("ObservingParameters/RApointing"));            
    originalDec = deg2rad(configParams.getDouble("ObservingParameters/DecPointing"));     
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
 * \brief Return the pointing coordinates (of the roll axis) of the spacecraft
 * 
 * \param time [s]
 *  
 * \return (alpha, delta)    RA & dec [rad] of the current pointing
 */

pair<double, double> Platform::getPointingCoordinates(double time)
{
    // Check if the request is going backwards in time. If so: complain.

    double timeInterval = time - internalTime;

    if (timeInterval < 0.0)
    {
        Log.warning("Platform: getPointingCoordinates() at time before previous request: Not Implemented. Returning current pointing coordinates.");
        return make_pair(currentRA, currentDec);
    }

    // Let the platfrom jitter until 'time'
    // Yaw, pitch, and roll are in [rad]

    double yaw, pitch, roll;
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

    return make_pair(currentRA, currentDec);
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
 * \param coordSC   (xSC, ySC, zSC): cartesian coordinates of the point in the spacecraft reference frame
 * \param useOriginalPointingCoordinates  If true: use original pointing coordinates (before jitter started)
 *                                        If false: use current pointing coordinates (affected by jitter)
 *
 * \return coordEQ  (xEQ, yEQ, zEQ): cartesian coordinates in the celestial equatorial reference frame
 */

arma::colvec Platform::spacecraftToEquatorialCoordinates(arma::colvec &coordSC, bool useOriginalPointingCoordinates)
{
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

    // Some handy abbreviations

    const double cosAlpha = cos(RA);
    const double sinAlpha = sin(RA);
    const double cosDelta = cos(dec);
    const double sinDelta = sin(dec);

    // The rotation matrices

    arma::mat R1 = {{cosAlpha, -sinAlpha, 0.0},
                    {sinAlpha,  cosAlpha, 0.0},
                    {     0.0,       0.0, 1.0}};

    arma::mat R2 = {{ sinDelta, 0.0, cosDelta},
                    {      0.0, 1.0,      0.0},
                    {-cosDelta, 0.0, sinDelta}};

    // The transformation

    arma::colvec coordEQ = R1 * R2 * coordSC;

    // That's it

    return coordEQ;
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

    arma::mat Ryaw = {{1.0,    0.0,     0.0},
                      {0.0, cosYaw, -sinYaw},
                      {0.0, sinYaw,  cosYaw}};

    arma::mat Rpitch = {{ cosPitch, 0.0, sinPitch},
                        {      0.0, 1.0,      0.0},
                        {-sinPitch, 0.0, cosPitch}};

    arma::mat Rroll = {{cosRoll, -sinRoll, 0.0},
                       {sinRoll,  cosRoll, 0.0},
                       {    0.0,      0.0, 1.0}}; 

    // Do the transformation

    return Ryaw  * Rpitch * Rroll * coord;
}



