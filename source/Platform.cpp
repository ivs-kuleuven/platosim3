
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
    originalRA  = deg2rad(configParams.getDouble("ObservingParameters/RApointing"));            
    originalDec = deg2rad(configParams.getDouble("ObservingParameters/DecPointing"));     
    currentRA   = originalRA;
    currentDec  = originalDec;
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
 * \brief Return the current sky pointing coordinates of the spacecraft
 * 
 * \return (alpha, delta)    RA & dec [rad] of the current pointing
 */

pair<double, double> Platform::getCurrentPointingCoordinates()
{
    return make_pair(currentRA, currentDec);
}









/**
 * \brief Update the pointing coordinates (varying due to jitter) until the given time point.
 * 
 * \param time [s]
 */

void Platform::updatePointingCoordinates(Telescope const &telescope, double time)
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
 *                 vector_after = R(yaw, pitch, roll) * vector_before
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



