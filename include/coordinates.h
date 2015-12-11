
#ifndef COORDINATES_H
#define COORDINATES_H

#include <tuple>
#include "units.h"


using namespace std;


// Strongly typed, but not scoped. Because CoordinateSystem::Equatorial is really long.

enum CoordinateSystem : short {Equatorial=0, Galactic=1, Ecliptic=2};



class Coordinates
{
    public:

        Coordinates(double longitude, double latitude, CoordinateSystem coordinateSystem=Equatorial);
        ~Coordinates(){}; 

        tuple<double, double> toEquatorial(Units units = Angle::degrees);
        tuple<double, double> toGalactic(Units units = Angle::degrees);
        tuple<double, double> toEcliptic(Units units = Angle::degrees);

 
        friend double angularDistanceBetween(Coordinates &coordinates1, Coordinates &coordinates2, Units units);


    protected:

        double longitude;                        // Equatorial, galactic, or ecliptic longitude [rad]
        double latitutde;                        // Equatorial, galactic, or ecliptic latitude  [rad]
 
    private:
    


};



#endif