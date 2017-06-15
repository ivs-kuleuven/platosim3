
#include "StarCatalog.h"
#include <valarray>




/**
 * \brief Default constructor
 */

 StarCatalog::StarCatalog()
 {

 }







/**
 * \brief Copy constructor
 */

StarCatalog::StarCatalog(const StarCatalog &starCatalog)
: starID(starCatalog.starID), RA(starCatalog.RA), dec(starCatalog.dec), Vmag(starCatalog.Vmag)
{

}









/**
 * \brief Destructor
 */

StarCatalog::~StarCatalog()
{

}







/**
 * \brief Return the number of stars in the catalog
 */


long StarCatalog::size()
{
    return starID.size();
}









/**
 * \brief Add the information of 1 more star in the catalog
 * 
 * \note There is no verification if the given starID is already in the catalog
 * 
 * \param starID     Star identification number      
 * \param RA         Right Ascension of the star [see angleUnit]
 * \param dec        Declination of the star [see angleUnit]
 * \param Vmag       Johnson V magnitude of the star
 * \param angleUnit  Angle::degrees if RA&dec are in degrees, Angle::radians if they are in radians.
 * 
 * \return None
 */

void StarCatalog::addStar(unsigned int starID, double RA, double dec, double Vmag, Unit angleUnit)
{
    this->starID.push_back(starID);
    this->RA.push_back(RA / angleUnit);
    this->dec.push_back(dec / angleUnit);
    this->Vmag.push_back(Vmag);
}










/**
 * \brief Operator [] for a StarCatalog. Returns a StarRecord of the 'index-th' star in the catalog.
 * 
 * \note  There is in-range checking for the index.
 * 
 * \param index  Integer between 0 and the number of stars in the catalog
 * \return       A StarRecord of the 'index-th' star in the catalog.
 */

StarRecord StarCatalog::operator[](unsigned int index) const
{
    // Check if the value of index is within the proper range

    if (index >= starID.size())
    {
        Log.error("StarCatalog[]: index out of range [0," + to_string(starID.size()) + "]");
        throw out_of_range("StarCatalog[]: index out of range");
    }

    // Return a StarRecord for this particular star

    return StarRecord(starID[index], RA[index], dec[index], Vmag[index]);
}










/**
 * \brief  Given a circle on the sky, return a new catalog with all stars from the current catalog within that circle.
 * 
 * \note   The stars right on the circle are also included in the catalog.
 * 
 * \param RA0              Right ascension of center point of the circle on the sky
 * \param dec0             Declination of center point of the circle on the sky
 * \param radius           Radius of the circle on the sky.
 * \param inputAngleUnit   Angle::degrees if input angles are in degrees, Angle::radians if in radians
 * 
 * \return            A StarCatalog with all the selected stars.
 */

StarCatalog StarCatalog::getStarsWithinRadiusFrom(double RA0, double dec0, double radius, Unit inputAngleUnit)
{
    // Create an empty star catalog

    StarCatalog newCatalog;

    // All computations are done in radians, so if RA0, dec0, and radius are expressed in degrees,
    // divide the degree unit away into radians.

    double RACircleCenter  = RA0    / inputAngleUnit;
    double decCircleCenter = dec0   / inputAngleUnit;
    double radiusCircle    = radius / inputAngleUnit;


    // Compute for each star in our catalog, the angular distance to (RA0, dec0)

    auto angularDistances = angularDistanceBetween(RACircleCenter, decCircleCenter, RA, dec, Angle::radians);

    // Copy those stars with a distance closer than 'radius' in the new catalog

    for (long n = 0; n < angularDistances.size(); ++n)
    {
        if (angularDistances[n] <= radiusCircle)
        {
            newCatalog.addStar(starID[n], RA[n], dec[n], Vmag[n], Angle::radians);
        }
    }

    return newCatalog;
}


/**
 * \brief  Calculate the apparent positions of the stars based on the current platform pointing coordinates.
 *
 * \detail
 *
 * This calculation is an approximation based on a circular earth orbit around the sun and *not* taking
 * the Lissajous orbit of the satellite around L2 into account. We do calculate the differential aberration
 * however which takes into account the aberration correction done for the Spacecraft pointing.
 * 
 * \param platform    the current platform from which the position of the Sun and the pointing coordinates are requested
 * 
 * \return            A StarCatalog with all the aberration corrected stars.
 */

StarCatalog StarCatalog::aberrate(Platform &platform, string aberrationCorrectionType, double startTime, double timeMiddle)
{
    using StringUtilities::dtos;

    // Create an empty star catalog

    StarCatalog newCatalog;

    //velocity direction of PLATO, assuming circular orbit in ecliptic plane with constant speed of 30 km/s, TODO: check the direction of rotation around the sun and adjust the sign of platoAngle accordingly
    double platoAngle = 2. * M_PI / 365. / 24. / 3600. * startTime;
    valarray<double> v = {cos(platoAngle), sin(platoAngle), 0.};

    //rotation matrix to compensate the aberration of light for the pointing direction, needed to calculate the differential aberration
    valarray<double> rot0 = {1., 0., 0.};
    valarray<double> rot1 = {0., 1., 0.};
    valarray<double> rot2 = {0., 0., 1.};

    //ratio of the velocity of PLATO to the speed of light
    constexpr double beta = 30. / 300000.;

    if (aberrationCorrectionType == "differential")
    {
        Log.info("StarCatalog::aberrate: applying differential aberration correction");

        // Request the current platform pointing coordinates (i.e. pointing of the Fast Camera's)

        double raPlatform, decPlatform;
        tie(raPlatform, decPlatform) = platform.getCurrentPointingCoordinates();

        double lambdaPlatform, betaPlatform;
        equatorial2ecliptic(raPlatform, decPlatform, lambdaPlatform, betaPlatform);

        //direction of the pointing
        valarray<double> p = {cos(lambdaPlatform) * cos(betaPlatform), sin(lambdaPlatform) * cos(betaPlatform), sin(betaPlatform)};

        //angle between velocity direction and pointing
        double pangle = acos((v * p).sum());

        //relativistically aberrated angle between velocity direction and pointing
        double oangle = atan2(sqrt(1. - beta * beta) * sin(pangle), cos(pangle) + beta);

        //rotation axis between velocity direction and pointing
        valarray<double> r = {p[1] * v[2] - p[2] * v[1], p[2] * v[0] - p[0] * v[2], p[0] * v[1] - p[1] * v[0]};
        r /= sqrt((r * r).sum()); 

        //rotation matrix for rotation axis r with angle difference after aberration, this reverses the aberration effect for the pointing direction
        double c = cos(oangle - pangle);
        double s = sin(oangle - pangle);
        double x = r[0], y = r[1], z = r[2];
        rot0 = {c + x * x * (1. - c), x * y * (1. - c) - z * s, x * z * (1. - c) + y * s};
        rot1 = {y * x * (1. - c) + z * s, c + y * y * (1. - c), y * z * (1. - c) - x * s};
        rot2 = {z * x * (1. - c) - y * s, z * y * (1. - c) + x * s, c + z * z * (1. - c)};

    }
    else
    {
        Log.info("StarCatalog::aberrate: applying absolute aberration correction");

    }

    for (long n = 0; n < starID.size(); ++n)
    {
        double raStar = RA[n];
        double decStar = dec[n];

        double lambdaStar, betaStar;
        equatorial2ecliptic(raStar, decStar, lambdaStar, betaStar);

        //direction of the star
        valarray<double> s = {cos(lambdaStar) * cos(betaStar), sin(lambdaStar) * cos(betaStar), sin(betaStar)};

        //angle between velocity direction and star direction
        double sangle = acos((v * s).sum());

        //relativistically aberrated angle between velocity direction and star direction
        double oangle = atan2(sqrt(1. - beta * beta) * sin(sangle), cos(sangle) + beta);

        //relativistically aberrated star direction
        valarray<double> a = s - v * cos(sangle);
        a = v * cos(oangle) + a / sqrt((a * a).sum()) * sin(oangle);

        //rotate aberrated star direction to compensate for aberrated pointing to get the differential aberrated star direction
        a = {(rot0 * a).sum(), (rot1 * a).sum(), (rot2 * a).sum()};

        //calculate ecliptic coordinates of aberrated star direction
        betaStar = atan(a[2] / sqrt(a[0] * a[0] + a[1] * a[1]));
        lambdaStar = atan2(a[1], a[0]);

        double raStarAberrated, decStarAberrated;
        ecliptic2equatorial(lambdaStar, betaStar, raStarAberrated, decStarAberrated);
        
        newCatalog.addStar(starID[n], raStarAberrated, decStarAberrated, Vmag[n], Angle::radians);

        // Write debugging info on the first star only

        if (n == 0)
        {
            Log.debug("StarCatalog::aberrate: ra[0], dec[0] = " + dtos(raStarAberrated, false, 8) + ", " + dtos(decStarAberrated, false, 8));
        }
    }

    return newCatalog;
}

