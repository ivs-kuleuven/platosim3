
#include "Camera.h"
#include "Units.h"
#include "Constants.h"


/**
 * \brief      Constructor
 *
 * \param configParam   Configuration parameters for the Camera
 * \param hdf5file      HFD5 file to write the camera information to
 * \param telescope     Telescope on which the camera is mounted on
 * \param sky           Sky object to query which stars we are seeing
 */
Camera::Camera(ConfigurationParameters &configParam, HDF5File &hdf5file, Telescope &telescope, Sky &sky)
: HDF5Writer(hdf5file), telescope(telescope), sky(sky), internalTime(0.0)
{
    // Parse the parameters from the configuration file.

    configure(configParam);
}





/**
 * @brief  Destructor
 */

Camera::~Camera()
{
}







/**
 * \brief      Expose the subField to the Sky, i.e. add flux to the detectors, 
 *             add Background and convolve with the PSF.
 *
 * \param[in]  detector  the Detector class
 */
void Camera::exposeDetector(Detector &detector)
{
    // auto starCatalog = sky.getStarsWithinRadiusFrom(alpha, delta, radius);  
    // double skyBackground = sky.getSkyBackground(alpha, delta)  

    // double tickInterval = telescope.getTickInterval();

    // while (currentTime < startingTime + exposureTime)
    // {
    //     telescope.updatePointingCoordinates(raOpticalAxis, decOpticalAxis, tickInterval);
    //     currentTime += tickInterval;

    //     for (auto star : starCatalog)
    //     {
    //         computeFocalPlaneCoordinates(star, Xmm, Ymm)
            
    //         if (subField.containsPoint(Xmm, Ymm))
    //         {
    //             subField.addFlux(Xmm, Ymm, flux);
    //         }
    //     }
    // }

    // subField.add(skyBackground * exposureTime);
    // subField.convolveWithPSF(psf);

    return;
}







/**
 * \brief Configure the Camera object using the ConfigurationParameters
 * 
 * \param configParam: the configuration parameters 
 */

void Camera::configure(ConfigurationParameters &configParam)
{
    plateScale            = configParam.getDouble("Camera/PlateScale");
    focalPlaneOrientation = deg2rad(configParam.getDouble("Camera/FocalPlaneOrientation"));
}







/**
 * \brief      select the PSF for the given star coordinates
 *
 * \param[in]  raStar   right ascension of the star [rad]
 * \param[in]  decStar  declination of the star     [rad]
 */

void Camera::selectPsf(double raStar, double decStar)
{
    pair<double, double> fpCoordinates = skyToFocalPlaneCoordinates(raStar, decStar);
}








/**
 * \brief Computes the (x,y) coordinates in the focal plane of a star 
 *        with given equatorial coordinates using a gnomonic projection
 *
 * \param raStar       Right ascension of the star [rad]
 * \param decStar      Declination of the star [rad]
 *
 * return pair (x,y):  Cartesian coordinate of the projected star in the focal plane in the FP-prime system [mm]
 */

pair<double, double> Camera::skyToFocalPlaneCoordinates(double raStar, double decStar)
{
    // Get the equatorial coordinates of the optical axis [rad]

    double raOpticalAxis, decOpticalAxis;
    tie(raOpticalAxis, decOpticalAxis) = telescope.getPointingCoordinates();

    // Project the sky to the focal plane in the "FP" coordinate system

    double denominator = cos(decOpticalAxis) * cos(decStar) * cos(raStar - raOpticalAxis) + sin(decOpticalAxis) * sin(decStar);
    double xFP = (sin(decOpticalAxis) * cos(decStar) * cos(raStar - raOpticalAxis) - cos(decOpticalAxis) * sin(decStar)) / denominator;
    double yFP =  cos(decStar) * sin(raStar - raOpticalAxis) / denominator;
    
    // Convert the FP coordinates into FP' coordinates 

    double xFPprime =  xFP * cos(focalPlaneOrientation) + yFP * sin(focalPlaneOrientation);
    double yFPprime = -xFP * sin(focalPlaneOrientation) + yFP * cos(focalPlaneOrientation);

    // Compute the conversion factor: conversion from [radians] to [mm].
    // Note that the plateScale is expressed in [arcsec/mm]

    double conversionFactor = 3600. * Constants::PI / 180.0 / plateScale;
    // double conversionFactor = 3600. * 180. / Constants::PI / 1000.0 / plateScale;

    Log.debug("Camera::skyToFocalPlaneCoordinates: conversionFactor = " + to_string(conversionFactor));

    // Return the scaled coordinates

    return make_pair(xFPprime * conversionFactor, yFPprime * conversionFactor);
}









/**
 * \brief Compute the equatorial sky coordinates of a star which has the given projected 
 *        focal plane (FP') coordinates (x,y)
 *        
 * \param xFPprimeStar     x-coordinate of the projected star in the focal plane in the FP-prime system [mm]
 * \param yFPprimeStar     y-coordinate of the projected star in the focal plane in the FP-prime system [mm]
 *
 * \return (alpha, delta)  Equatorial coordinates (RA & Dec) of the star [rad]
 */

pair<double, double> Camera::focalPlaneToSkyCoordinates(double xFPprimeStar, double yFPprimeStar)
{    
    if (xFPprimeStar == 0.0 and yFPprimeStar == 0.0)
    {
        return make_pair(0.0, 0.0);
    }
    
    // Compute the conversion factor from [mm/radian] to [arcsec/pixel] and convert the
    // focal plane coordinates from [mm] to [rad].

    double conversionFactor = 3600. * Constants::PI / 180.0 / plateScale;
    double xFPprime = xFPprimeStar / conversionFactor;
    double yFPprime = yFPprimeStar / conversionFactor;

    // Convert the FP' coordinates into FP coordinates

    double xFP =  xFPprime * cos(focalPlaneOrientation) - yFPprime * sin(focalPlaneOrientation);
    double yFP =  xFPprime * sin(focalPlaneOrientation) + yFPprime * cos(focalPlaneOrientation);

    // Get the equatorial coordinates of the optical axis [rad]

    double raOpticalAxis, decOpticalAxis;
    tie(raOpticalAxis, decOpticalAxis) = telescope.getPointingCoordinates();

    // Project the focal plane in the "FP" coordinate system to the sky

    double rho = sqrt(xFP*xFP+yFP*yFP);
    double c = atan(rho);
    double decStar = asin(cos(c)*sin(decOpticalAxis)+(-xFP*sin(c)*cos(decOpticalAxis))/rho);
    double raStar = raOpticalAxis + atan2(yFP*sin(c), rho*cos(decOpticalAxis)*cos(c)+xFP*sin(decOpticalAxis)*sin(c));
    

    // Return the equatorial coordinates

    return make_pair(raStar, decStar);
}
