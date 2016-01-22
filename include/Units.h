
#ifndef UNITS_H
#define UNITS_H


typedef double Unit;


struct Angle
{
    static constexpr double radians = 1.0;                             // default unit
    static constexpr double degrees = 57.29577951308232;               // 180/pi
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


#endif
