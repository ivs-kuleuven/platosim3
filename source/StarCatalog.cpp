
#include "StarCatalog.h"




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
 * \brief Move constructor
 */

StarCatalog::StarCatalog(StarCatalog &&starCatalog)
: starID(move(starCatalog.starID)), RA(move(starCatalog.RA)), dec(move(starCatalog.dec)), Vmag(move(starCatalog.Vmag))
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


