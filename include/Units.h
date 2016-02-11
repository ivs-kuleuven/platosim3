
#ifndef UNITS_H
#define UNITS_H


typedef double Unit;


struct Angle
{
    static constexpr double radians = 1.0;                             // default unit
    static constexpr double degrees = 57.29577951308232;               // 180/pi
};



struct SolidAngle
{
    static constexpr double steradians = 1.0;                          // default unit
    static constexpr double squareDegrees = 3282.8063500117437;        // (180/pi)^2
};





/**
 * @brief convert angle units from degrees to radians
 * 
 * @param angleInDegrees: angle in degrees
 * @return angle in radians
 */

inline double deg2rad(const double angleInDegrees)
{
    return angleInDegrees / Angle::degrees;
}







/**
 * @brief convert angle units from radians to degrees
 * 
 * @param angleInRadians: angle in radians
 * @return: angle in degrees
 */

inline double rad2deg(const double angleInRadians)
{
    return angleInRadians * Angle::degrees;
}






/**
 * @brief  Convert solid angle units from square degrees to steradians
 * 
 * @param solidAngleInSquareDegrees  Solid angle in square degrees
 * @return  Solid angle in steradians
 */

inline double sqDeg2sr(const double solidAngleInSquareDegrees)
{
    return solidAngleInSquareDegrees / SolidAngle::squareDegrees;
}







/**
 * @brief  Convert solid angle units from steradians to square degrees
 * 
 * @param angleInSteradians   Solid angle in steradians
 * @return                    Solid angle in square degrees
 */

inline double sr2sqDeg(const double angleInSteradians)
{
    return angleInSteradians * SolidAngle::squareDegrees;
}








#endif
