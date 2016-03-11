#include "gtest/gtest.h"

#include "gtest_definitions.h"

#include "Units.h"
#include "Constants.h"
#include "Logger.h"
#include "ArrayOperations.h"

using namespace std;


void checkArraysToBeEqual(arma::fmat arr1, arma::fmat arr2);





// TODO: * Add tests for other rotation angles
//       * Add a test where we create a round image with some defined data and every pixel around that is black.
//         Then rotate that image by different degrees and verify that some parameters are unchanged, e.g.
//         the average pixel values, the sum of all the pixel values (flux conservative) etc.
//       * Create a series of rotations which add up to 360 degrees and then compare to the original.
//       * Rotate by 90 degrees and -270 degrees and compare the results, use different complementary angles to do the same.
//       
TEST(ArrayOperationsTest, Rotation)
{

    LOG_STARTING_OF_TEST

    arma::fmat bigNulledArray = {
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

    // Test array where there are no zero (0) values at the sides.
    arma::fmat smallArray = {
        {1., 2., 3.},
        {4., 5., 6.},
        {7., 8., 9.}
    };

    // Expected result after rotating the smallArray by 90 degrees
    arma::fmat rotatedSmallArray90 = {
        {3., 6., 9.},
        {2., 5., 8.},
        {1., 4., 7.}
    };

    // Test array where zeros (0) are added
    arma::fmat smallNulledArray = {
        {0. , 0., 0., 0., 0.},
        {0. , 1., 2., 3., 0.},
        {0. , 4., 5., 6., 0.},
        {0. , 7., 8., 9., 0.},
        {0. , 0., 0., 0., 0.}
    };

    // Expected result after rotating the smallNulledArray by 90 degrees
    arma::fmat rotatedSmallNulledArray90 = {
        {0. , 0., 0., 0., 0.},
        {0. , 3., 6., 9., 0.},
        {0. , 2., 5., 8., 0.},
        {0. , 1., 4., 7., 0.},
        {0. , 0., 0., 0., 0.}
    };


    
    arma::fmat arr = rotateArray(smallArray, deg2rad(90.0));
    checkArraysToBeEqual(rotatedSmallArray90, arr);

//    printArray(smallArray, "Original smallArray");
//    printArray(arr, "Rotated smallArray");
    




    arr = rotateArray(smallNulledArray, deg2rad(90.0));
    checkArraysToBeEqual(rotatedSmallNulledArray90, arr);

//    printArray(smallNulledArray, "Original smallNulledArray");
//    printArray(arr, "Rotated smallNulledArray");




    arr = rotateArray(bigNulledArray, deg2rad(90.0));
    EXPECT_FLOAT_EQ(50.0, arr(9, 9));
    EXPECT_FLOAT_EQ(10.0, arr(5, 5));
    EXPECT_FLOAT_EQ(10.0, arr(13, 13));
    EXPECT_FLOAT_EQ(20.0, arr(13, 5));
    EXPECT_FLOAT_EQ(40.0, arr(5, 13));


//    printArray(bigNulledArray, "Original bigNulledArray");
//    printArray(arr, "Rotated bigNulledArray");

}



void checkArraysToBeEqual(arma::fmat arr1, arma::fmat arr2)
{
    EXPECT_EQ(arr1.n_rows, arr2.n_rows);
    EXPECT_EQ(arr1.n_cols, arr2.n_cols);

    for (unsigned int ix=0; ix<arr1.n_rows; ix++)
    {
        for (unsigned int jy=0; jy<arr1.n_cols; jy++)
        {
            EXPECT_NEAR(arr1(ix, jy), arr2(ix, jy), 0.0001);
        }
    }
}