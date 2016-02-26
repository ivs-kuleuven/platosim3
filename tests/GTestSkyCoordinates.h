#include "gtest/gtest.h"

#include "gtest_definitions.h"

#include "Units.h"
#include "Constants.h"
#include "SkyCoordinates.h"

using namespace std;

// TODO: Find some good boundary tests for calculating the angular distance bewteen two coordinates


TEST(SkyCoordinatesTest, angularDistance)
{
    LOG_STARTING_OF_TEST

    // Test the angular distance between two points on the sky.
    // Comparison values were computed using http://www.asdc.asi.it/dist.html

    SkyCoordinates star1(0.0, 0.0, Angle::degrees);
    SkyCoordinates star2(178, -70, Angle::degrees);

    double angle = angularDistanceBetween(star1, star2, Angle::degrees);
    ASSERT_NEAR(angle, 109.98730, 1.e-5);

}







TEST(SkyCoordinates, conversionEquatorialToGalactic)
{
    LOG_STARTING_OF_TEST

    double longitude, latitude;

    // Test the conversion for different coordinates in degrees
    // Comparison values were computed using https://ned.ipac.caltech.edu/forms/calculator.html
    // Note: Epoch=1950.0. Website's RA is in h:m:s.

    SkyCoordinates star1(0.0, 0.0,  Angle::degrees);
    tie(longitude, latitude) = star1.toGalactic(Angle::degrees);
    ASSERT_NEAR(longitude, 97.74216087, 1.e-4);
    ASSERT_NEAR(latitude, -60.18102400, 1.e-4);

    SkyCoordinates star2(120.0, -45.0,  Angle::degrees);
    tie(longitude, latitude) = star2.toGalactic(Angle::degrees);
    ASSERT_NEAR(longitude, 260.18944349, 1.e-4);
    ASSERT_NEAR(latitude,   -7.70210503, 1.e-4);

    SkyCoordinates star3(180.0,  90.0,  Angle::degrees);
    tie(longitude, latitude) = star3.toGalactic(Angle::degrees);
    ASSERT_NEAR(longitude, 123.00000000, 1.e-4);
    ASSERT_NEAR(latitude,   27.40000000, 1.e-4);


    // Do the same testing but in radians rather than degrees

    SkyCoordinates star4(deg2rad(0.0), deg2rad(0.0),  Angle::radians);
    tie(longitude, latitude) = star4.toGalactic(Angle::radians);
    ASSERT_NEAR(longitude, deg2rad(97.74216087),  deg2rad(1.e-4));
    ASSERT_NEAR(latitude,  deg2rad(-60.18102400), deg2rad(1.e-4));

    SkyCoordinates star5(deg2rad(120.0), deg2rad(-45.0),  Angle::radians);
    tie(longitude, latitude) = star5.toGalactic(Angle::radians);
    ASSERT_NEAR(longitude, deg2rad(260.18944349), deg2rad(1.e-4));
    ASSERT_NEAR(latitude,   deg2rad(-7.70210503), deg2rad(1.e-4));

    SkyCoordinates star6(deg2rad(180.0),  deg2rad(90.0),  Angle::radians);
    tie(longitude, latitude) = star6.toGalactic(Angle::radians);
    ASSERT_NEAR(longitude, deg2rad(123.00000000), deg2rad(1.e-4));
    ASSERT_NEAR(latitude,   deg2rad(27.40000000), deg2rad(1.e-4));

       
}





TEST(SkyCoordinates, conversionEquatorialToEcliptic)
{
    LOG_STARTING_OF_TEST

    double longitude, latitude;

    // Test the conversion for different coordinates in degrees
    // Comparison values were computed using https://ned.ipac.caltech.edu/forms/calculator.html
    // Note: Epoch=2000.0. Website's RA is in h:m:s.

    SkyCoordinates star1(0.0, 0.0,  Angle::degrees);
    tie(longitude, latitude) = star1.toEcliptic(Angle::degrees);
    ASSERT_NEAR(longitude, 0.0, 1.e-4);
    ASSERT_NEAR(latitude,  0.0, 1.e-4);

    SkyCoordinates star2(120.0, -45.0,  Angle::degrees);
    tie(longitude, latitude) = star2.toEcliptic(Angle::degrees);
    ASSERT_NEAR(longitude, 141.56549538, 5.e-4);
    ASSERT_NEAR(latitude,  -63.16948415, 5.e-4);

    SkyCoordinates star3(180.0,  90.0,  Angle::degrees);
    tie(longitude, latitude) = star3.toEcliptic(Angle::degrees);
    ASSERT_NEAR(longitude, 90.00000, 5.e-4);
    ASSERT_NEAR(latitude,  66.56071, 5.e-4);


    // Do the same testing but in radians rather than degrees

    SkyCoordinates star4(deg2rad(0.0), deg2rad(0.0),  Angle::radians);
    tie(longitude, latitude) = star4.toEcliptic(Angle::radians);
    ASSERT_NEAR(longitude, deg2rad(0.0), deg2rad(1.e-4));
    ASSERT_NEAR(latitude,  deg2rad(0.0), deg2rad(1.e-4));

    SkyCoordinates star5(deg2rad(120.0), deg2rad(-45.0),  Angle::radians);
    tie(longitude, latitude) = star5.toEcliptic(Angle::radians);
    ASSERT_NEAR(longitude, deg2rad(141.56549538), deg2rad(5.e-4));
    ASSERT_NEAR(latitude,  deg2rad(-63.16948415), deg2rad(5.e-4));

    SkyCoordinates star6(deg2rad(180.0),  deg2rad(90.0),  Angle::radians);
    tie(longitude, latitude) = star6.toEcliptic(Angle::radians);
    ASSERT_NEAR(longitude, deg2rad(90.00000), deg2rad(5.e-4));
    ASSERT_NEAR(latitude,  deg2rad(66.56071), deg2rad(5.e-4)); 

}