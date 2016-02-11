
#include "Platform.h"



/**
 * \brief Constructor
 * 
 * \param configurationParameters  Configuration parameters as read from the (e.g. yaml) inputfile
 * \param hdf5File                 HDF5 file where to write any output
 * \param jitterGenerator          Generator of the yaw, pitch roll variations due to spacecraft jitter
 */

Platform::Platform(ConfigurationParameters configParams, HDF5File &hdf5File, JitterGenerator &jitterGenerator)
: HDF5Writer(hdf5File), jitterGenerator(jitterGenerator)
{
    // Configure the Platfrom object

    configure(configParams);
}







/**
 * \brief Destructor
 */

Platform::~Platform()
{

}








/**
 * \brief Configure the Platform object
 * 
 * \param configParams  The configuration parameters from the input parameters file
 */

void Platform::configure(ConfigurationParameters &configParams)
{
    currentRA  = deg2rad(configParams.getDouble("ObservingParameters/RApointing"));            
    currentDec = deg2rad(configParams.getDouble("ObservingParameters/DecPointing"));     
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










/**
 * \brief Transforms the coordinates of a vector after yaw, pitch, and roll rotations.
 * 
 * \details We define the yaw, the pitch, and the roll as rotation angles around respectively 
 *          the X_SC, Y_SC, and Z_SC axes, (SC standing for spacecraft reference frame) such that 
 *          the angles increase with a clockwise rotation, when looking along the positive axes. 
 *          First a roll rotation is done around the Z_SC axis, then a pitch rotation is done 
 *          around the rotated Y_SC axis, and finally a yaw rotation is done around the already
 *          twice-rotated X_SC axis. The combined rotation matrix is: 
 *                 R(yaw, pitch, roll) = R(yaw) * R(pitch) * R(roll)
 *          Which transforms
 *                 vector_before = R(yaw, pitch, roll) * vector_after
 *          where vector_before and vector_after are the coordinates of the same point before and
 *          after the yaw.pitch,roll rotations. The inverse transformation (usd in this function) is
 *                 R^{-1}(yaw, pitch, roll) = R^t(roll) * R^t(pitch) * R^t(yaw)
 *          which transforms
 *                 vector_after = R(yaw, pitch, roll) * vector_before
 *          Here ^t denotes the transpose. 
 * 
 * \param coord Coordinates of the point in the spacecraft reference frame before the rotation 
 * \param yaw   Yaw angle   [rad]
 * \param pitch Pitch angle [rad]
 * \param roll  Roll angle  [rad]
 * \return      Coordinates of the point in the spacecraft reference frame after the rotations
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

    return Rroll.t() * Rpitch.t() * Ryaw.t() * coord;
}



