/**
 * \class     Camera 
 * 
 * \brief     Handle the distortions and effects that are cause by the optical system of the Telescope.
 * 
 * \details
 * 
 * The Camera is basically the set of six lenses with their mechanical mounts and support 
 * structure also known as the TOU or Telescope Optical Units. The lenses distort the 
 * incoming light in several ways. The following effects are due to the setup and 
 * characteristics of the lenses:
 *
 *   * Image Quality (Enclosed Energy)
 *   * Optical distortion
 *   * Vignetting
 *   * Point Spread Function (PSF)
 *   * PSF Breathing due to thermal variations
 *   * Transmission Efficiency
 *   * Straylight
 *   * Lens degradation and contamination (??)
 * 
 * Not all above effects are implemented in the PLATO Simulator at this point. We concentrate on the 
 * most distinct effects like PSF, optical distortion, and vignetting.
 * 
 * The lenses are the main source of point source spreading over the detector array. The camera is 
 * therefore the obvious choice for applying the PSF correction. The PSF itself is described in it’s 
 * own class, see PointSpreadFunction Class.
 * 
 * 
 * 
 */
#include "Camera.h"
#include "Units.h"
#include "Constants.h"
#include "PointSpreadFunction.h"

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

    // Initialize and load the PSF. This will open the PSF HDF5 file and perform some basic checking, 
    // but the psf itself will only be loaded by the psf.select(radius) method.

    psf = new PointSpreadFunction(configParam);
}






/**
 * \brief  Destructor
 */

Camera::~Camera()
{
    delete psf;
}








/**
 * \brief      Expose the subField to the Sky, i.e. add flux to the detectors,
 *             add Background and convolve with the PSF.
 *
 * \param[in]  detector      the Detector class
 * \param[in]  startTime     start time of the exposure [seconds]
 * \param[in]  exposureTime  duration of one exposure [seconds]
 */
void Camera::exposeDetector(Detector &detector, double startTime, double exposureTime)
{
    // Get the focal plane coordinates of the center of the subfield (in [mm]), 
    // and the diagonal length of the subfield (converted from [mm] to [rad]).
    // These quantities are fixed, i.e. independent of any jitter.
    // Note: diagonalLength is in [mm], platescale is in [arcsec/mm], radius is in [rad]

    double Xmm, Ymm;
    tie(Xmm, Ymm) = detector.getFocalPlaneCoordinatesOfSubfieldCenter();
    double diagonalLength = detector.getDiagonalLengthOfSubfield();
    double radius = deg2rad(diagonalLength / 2.0 * plateScale / 3600.);

    // Compute the (alpha, delta) equatorial coordinates in [rad] of the center of the subfield [rad]

    double RA, dec;
    tie(RA, dec) = focalPlaneToSkyCoordinates(Xmm, Ymm);

    // Get a catalog of stars that fall on the subfield. Take the radius a bit larger so that the 
    // queried area includes possible small shifts of the projected subfield because of jitter.

    auto starCatalog = sky.getStarsWithinRadiusFrom(RA, dec, radius * 1.1, Angle::radians);  

    // If the telescope and/or platform show small variations (e.g. due to jitter) during the exposure,
    // the exposure time is split up in many small intervals, to track the effect of these variations
    // on the exposure. The largest time interval for which the variations can still be reliably tracked
    // is called the heartbeatInterval. The time step used should be either the heartbeat interval or the
    // expsosure time whatever is smallest. 

    double timeStep = min(telescope.getHeartbeatInterval(), exposureTime);

    // Later on we will have to convert from magnitudes to fluxes. Precompute a constant prefactor.
    // 100023.8 is the photon flux [photons/s/cm^2/nm] for a V=0 G2V-star.
    // Units of fluxFactor: [photons/s]
  
    const double fluxFactor = 10023.8 * throughputBandwidth * telescope.getTransmissionEfficiency() * telescope.getLightCollectingArea(); 

    // Update the internal clock

    internalTime = startTime;

    // Take the flux of point sources (stars) into account.
    // Break up the exposure time in small intervals (hearbeat intervals) to track jitter while exposing.

    while (internalTime < startTime + exposureTime)
    {
        // Update the clock. Normally with 'timeStep', but if adding timeStep would overstep
        // the total exposure time, take the small rest time instead.

        timeStep = min(timeStep, startTime + exposureTime - internalTime);
        internalTime += timeStep;

        // Let the telescope pointing evolve over a small time interval

        telescope.updatePointingCoordinates(internalTime);

        // Loop over all stars in the catalog, and add their flux to the subfield

        for (int n = 0; n < starCatalog.size(); n++)
        {
            // Get the focal plane coordinates (in [mm]) of this particular star
            
            auto star = starCatalog[n];
            tie(Xmm, Ymm) = skyToFocalPlaneCoordinates(star.RA, star.dec);

            // Compute the flux [photons] of this star

            double flux = fluxFactor * pow(10.0, -0.4 * star.Vmag) * timeStep;

            // Add the flux to the detector. The latter checks if the star falls on a pixel of the subfield.

            detector.addFlux(Xmm, Ymm, flux);
        }
    }

    // Take the flux of the stellar background and the zodiacal light into account.
    // Use one value for the entire subfield.

    const double energyOfOnePhoton = Constants::CLIGHT * Constants::HPLANCK / (throughputLambdaC * 1.e-9);                       // [J]
    const double lambda1 = (throughputLambdaC - throughputBandwidth/2.0) * 1.e-9;                                                // [m]
    const double lambda2 = (throughputLambdaC + throughputBandwidth/2.0) * 1.e-9;                                                // [m]
    
    double backgroundFlux = (sky.zodiacalFlux(RA, dec, lambda1, lambda2) + sky.stellarBackgroundFlux(RA, dec, lambda1, lambda2))
                            * exposureTime * telescope.getTransmissionEfficiency() * telescope.getLightCollectingArea()
                            * telescope.getFOVsolidAngle() / energyOfOnePhoton;                                                  // [photons/exposure]
    
    detector.addFlux(backgroundFlux);


    // Convolve with the point spread function

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
    throughputBandwidth   = configParam.getDouble("Camera/ThroughputBandwidth");
    throughputLambdaC     = configParam.getDouble("Camera/ThroughputLambdaC");
}








/**
 * \brief      select the PSF for the given star coordinates.
 * 
 * \details
 * 
 * This method selects and rotates the PSF.
 *
 * \param[in]  raStar   right ascension of the star [rad]
 * \param[in]  decStar  declination of the star     [rad]
 */

void Camera::selectPsf(double raStar, double decStar)
{

    // Get the equatorial coordinates of the optical axis [rad]

    double raOpticalAxis, decOpticalAxis;
    tie(raOpticalAxis, decOpticalAxis) = telescope.getCurrentPointingCoordinates();

    // Calculate the angular separation (in [radians]) between the star and the optical axis.
    // Use that angle to select the proper PSF.

    Coordinates opticalAxis(raOpticalAxis, decOpticalAxis, Angle::radians);
    Coordinates star(raStar, decStar, Angle::radians);

    double radius = angularDistanceBetween(opticalAxis, star, Angle::radians);

    psf->select(radius);

    // Calculate the rotation angle

    double xFP, yFP;
    tie(xFP, yFP) = skyToNormalizedFocalPlaneCoordinates(raStar, decStar);

    double angle = atan2(yFP, xFP);

    psf->rotate(angle);

}














/**
 * \brief      Calculate the gnomonic radial distance with respect to the optical axis in the focal plane
 *
 * \param[in]  xFPprime Cartesian x-coordinate of the projected star in the focal plane in the FP-prime system [mm]
 * \param[in]  yFPprime Cartesian y-coordinate of the projected star in the focal plane in the FP-prime system [mm]
 *
 * \return     the field radial distance (gnomonic) with respect to the line of sight in the sky [deg]
 */
double Camera::getGnomonicRadialDistanceFromOpticalAxis(double xFPprime, double yFPprime)
{
    // Convert from planar focal plane coordinates [mm] to angular focal plane coordinates [radians]
    // We don't need the angular FP coordinates themselves, just the tangens of them. 
    // So we avoid taking first the arctan() and then the tan()

    const double focalLength = 3600. * 180.0 / Constants::PI / 1000.0 / plateScale;
    const double tanx = xFPprime / focalLength;               // = tan(atan(xFPprime / focalLength))
    const double tany = yFPprime / focalLength;

    return rad2deg(acos(1.0/sqrt(1.0 + tanx*tanx + tany*tany)));
}




/**
 * @brief      Calculate the gnomonic radial distance with respect to the optical axis in the normalized focal plane
 *
 * @param[in]  xFPprime  normalized focal plane x-coordinate [rad]
 * @param[in]  yFPprime  normalized focal plane y-coordinate [rad]
 *
 * @return     the field radial distance (gnomonic) with respect to the line of sight in the sky [rad]
 */
double Camera::getGnomonicRadialDistanceFromOpticalAxisNormalized(double xFPprime, double yFPprime)
{
    double tanx = tan(xFPprime);
    double tany = tan(yFPprime);

    return acos(1.0/sqrt(1.0 + tanx*tanx + tany*tany));
}





/**
 * \brief Computes the (x,y) coordinates in the focal plane of a star with given equatorial coordinates
 *        using a gnomonic projection.
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
    tie(raOpticalAxis, decOpticalAxis) = telescope.getCurrentPointingCoordinates();

    // Project the sky to the focal plane in the "FP" coordinate system (gnomonic projection)

    const double denominator = cos(decOpticalAxis) * cos(decStar) * cos(raStar - raOpticalAxis) + sin(decOpticalAxis) * sin(decStar);
    const double xFPrad= ( sin(decOpticalAxis) * cos(decStar) * cos(raStar - raOpticalAxis) - cos(decOpticalAxis) * sin(decStar)) / denominator;
    const double yFPrad =  cos(decStar) * sin(raStar - raOpticalAxis) / denominator;

    // Convert the FP coordinates into FP' coordinates. Both are in [rad]. 

    const double xFPprime =  xFPrad * cos(focalPlaneOrientation) + yFPrad * sin(focalPlaneOrientation);
    const double yFPprime = -xFPrad * sin(focalPlaneOrientation) + yFPrad * cos(focalPlaneOrientation);

    // Convert from angular focal plane coordinates [radians] to planar focal plane coordinates [mm].
    // Note that the plateScale is expressed in [arcsec/um]

    const double focalLength = 3600. * 180.0 / Constants::PI / 1000.0 / plateScale;
    const double xFPprime_mm = focalLength * tan(xFPprime);
    const double yFPprime_mm = focalLength * tan(yFPprime);

    // Return the scaled coordinates

    return make_pair(xFPprime_mm, yFPprime_mm);
}






/**
 * \brief Computes the (x,y) coordinates in the normalized focal plane of a star with given equatorial coordinates
 *        using a gnomonic projection.
 *
 * \param raStar       Right ascension of the star [rad]
 * \param decStar      Declination of the star [rad]
 *
 * return pair (x,y):  Cartesian coordinate of the projected star in the normalized focal plane in the FP-prime system [radians]
 */

pair<double, double> Camera::skyToNormalizedFocalPlaneCoordinates(double raStar, double decStar)
{
    // Get the equatorial coordinates of the optical axis [rad]

    double raOpticalAxis, decOpticalAxis;
    tie(raOpticalAxis, decOpticalAxis) = telescope.getCurrentPointingCoordinates();

    // Project the sky to the focal plane in the "FP" coordinate system (gnomonic projection)

    double denominator = cos(decOpticalAxis) * cos(decStar) * cos(raStar - raOpticalAxis) + sin(decOpticalAxis) * sin(decStar);
    double xFP = ( sin(decOpticalAxis) * cos(decStar) * cos(raStar - raOpticalAxis) - cos(decOpticalAxis) * sin(decStar)) / denominator;
    double yFP =  cos(decStar) * sin(raStar - raOpticalAxis) / denominator;

    // Convert the FP coordinates into FP' coordinates 

    double xFPprime =  xFP * cos(focalPlaneOrientation) + yFP * sin(focalPlaneOrientation);
    double yFPprime = -xFP * sin(focalPlaneOrientation) + yFP * cos(focalPlaneOrientation);

    // Return the scaled coordinates

    return make_pair(xFPprime, yFPprime);
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

    // Convert from planar focal plane coordinates [mm] to angular focal plane coordinates [radians] 

    const double focalLength = 3600. * 180.0 / Constants::PI / 1000.0 / plateScale;
    const double xFPprime = atan(xFPprimeStar / focalLength);
    const double yFPprime = atan(yFPprimeStar / focalLength);

    // Convert from the FP' to the FP reference system.

    const double xFP =  xFPprime * cos(focalPlaneOrientation) - yFPprime * sin(focalPlaneOrientation);
    const double yFP =  xFPprime * sin(focalPlaneOrientation) + yFPprime * cos(focalPlaneOrientation);


    // Get the equatorial coordinates of the optical axis [rad]

    double raOpticalAxis, decOpticalAxis;
    tie(raOpticalAxis, decOpticalAxis) = telescope.getCurrentPointingCoordinates();

    // Project the focal plane in the "FP" coordinate system to the sky

    const double rho = sqrt(xFP*xFP+yFP*yFP);
    const double c = atan(rho);
    const double raStar = raOpticalAxis + atan2(yFP * sin(c), rho * cos(decOpticalAxis) * cos(c) + xFP * sin(decOpticalAxis) * sin(c));
    const double decStar = asin(cos(c) * sin(decOpticalAxis) - (xFP * sin(c) * cos(decOpticalAxis)) / rho);


    // Return the equatorial coordinates

    return make_pair(raStar, decStar);
}
