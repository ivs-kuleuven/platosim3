#include "gtest/gtest.h"

#include "PointSpreadFunction.h"
#include "Units.h"

using namespace std;

TEST(PointSpreadFunctionTest, Constructor_ConfigurationParameters)
{
    ConfigurationParameters cp = ConfigurationParameters("../testData/input_PointSpreadFunctionTest.yaml");

    PointSpreadFunction psf = PointSpreadFunction(cp);

}

TEST(PointSpreadFunctionTest, DISABLED_Rotation)
{
    unsigned int numRowsPsfMap = 19;
    unsigned int numColumnsPsfMap = 19;
    double angle = 90.0;

    arma::Mat<float> psfMap = {
               {0., 0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0., 0., 0.}, \
               {0., 0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0., 0., 0.}, \
               {0., 0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0., 0., 0.}, \
               {0., 0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0., 0., 0.}, \
               {0., 0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0., 0., 0.}, \
               {0., 0., 0., 0., 0., 20., 1., 1., 1.,  1., 1., 1., 1., 10., 0., 0., 0., 0., 0.}, \
               {0., 0., 0., 0., 0.,  1., 1., 1., 1.,  1., 1., 1., 1.,  1., 0., 0., 0., 0., 0.}, \
               {0., 0., 0., 0., 0.,  1., 1., 1., 1.,  1., 1., 1., 1.,  1., 0., 0., 0., 0., 0.}, \
               {0., 0., 0., 0., 0.,  1., 1., 1., 1.,  1., 1., 1., 1.,  1., 0., 0., 0., 0., 0.}, \
               {0., 0., 0., 0., 0.,  1., 1., 1., 1., 50., 1., 1., 1.,  1., 0., 0., 0., 0., 0.}, \
               {0., 0., 0., 0., 0.,  1., 1., 1., 1.,  1., 1., 1., 1.,  1., 0., 0., 0., 0., 0.}, \
               {0., 0., 0., 0., 0.,  1., 1., 1., 1.,  1., 1., 1., 1.,  1., 0., 0., 0., 0., 0.}, \
               {0., 0., 0., 0., 0.,  1., 1., 1., 1.,  1., 1., 1., 1.,  1., 0., 0., 0., 0., 0.}, \
               {0., 0., 0., 0., 0., 10., 1., 1., 1.,  1., 1., 1., 1., 40., 0., 0., 0., 0., 0.}, \
               {0., 0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0., 0., 0.}, \
               {0., 0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0., 0., 0.}, \
               {0., 0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0., 0., 0.}, \
               {0., 0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0., 0., 0.}, \
               {0., 0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0.,  0., 0., 0., 0., 0., 0.}, \
    };

    arma::Mat<float> rotatedPsfMap(arma::size(psfMap), arma::fill::zeros);

    float background = 0.0;  // this is the background color
    float rads = deg2rad(angle);
    float cs = cos(-rads);  // precalculate these values
    float ss = sin(-rads);
    float xcenter = (float)(numRowsPsfMap)/2.0;   // use float here!
    float ycenter = (float)(numColumnsPsfMap)/2.0;

    for (int row = 0; row < numRowsPsfMap; ++row) {
       for (int column = 0; column < numColumnsPsfMap; ++column) {
          // now find the pixel of the original image that is rotated to (row, column)
          // rotation formula assumes that origin = top-left and y points down
          int rorig = ycenter + ((float)(row)-ycenter)*cs - ((float)(column)-xcenter)*ss;
          int corig = xcenter + ((float)(row)-ycenter)*ss + ((float)(column)-xcenter)*cs;
          // now get the pixel value if you can
          float pixel = background; // in case there is no original pixel
          if (rorig >= 0 && rorig < numRowsPsfMap && corig >= 0 && corig < numColumnsPsfMap) {
             pixel = psfMap(rorig, corig);
          }
          rotatedPsfMap(row, column) = pixel;
       }
    }

}