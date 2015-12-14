
#ifndef UNITS_H
#define UNITS_H


typedef double Units;


struct Angle
{
    static constexpr double radians = 1.0;                             // default unit
    static constexpr double degrees = 57.29577951308232;               // 180/pi
};


#endif