
#ifndef COORDINATES_H
#define COORDINATES_H

#include <cmath>

#include "Units.h"
#include "Logger.h"


using namespace std;




class Coordinates
{
    public:

        Coordinates(double RA, double decl, Unit angleUnit = Angle::degrees);
        ~Coordinates(); 

        pair<double, double> toGalactic(Unit angleUnit = Angle::degrees);
        pair<double, double> toEcliptic(Unit angleUnit = Angle::degrees);

 
        friend double angularDistanceBetween(Coordinates &coordinates1, Coordinates &coordinates2, Unit angleUnit);
        friend vector<double> angularDistanceBetween(const double RA0, const double dec0, const vector<double> &RA, const vector<double> &dec, Unit angleUnit);

    protected:

        void equatorial2ecliptic(const double alpha, const double delta, double &lambda, double &beta);
        void equatorial2galactic(const double alpha, const double delta, double &l, double &b);

        double RA;                    // Equatorial right ascension [rad]
        double decl;                  // Equatorial declination     [rad]
 
    private:

        const double obliquity;              // Obliquity of the ecliptic = 23.439 deg  [rad]
        const double inclGalPlane;           // Inclination of the galactic plane = 62.6 deg in B1950 [rad]
        const double alphaN;                 // Right ascension of the ascending node of the galactic plane = 282.25 deg in B1950 [rad]
        const double l0;                     // Galactic longitude of the ascending node of the galactic plane = 33 deg in B1950 [rad]
    
};



#endif