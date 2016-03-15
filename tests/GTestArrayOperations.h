#include "gtest/gtest.h"

#include "gtest_definitions.h"

#include "Units.h"
#include "Constants.h"
#include "Logger.h"
#include "ArrayOperations.h"
#include "StringUtilities.h"

using namespace std;
using namespace ArrayOperations;


void checkArraysToBeEqual(arma::fmat arr1, arma::fmat arr2);






// This test is actually to guide and practice the access of elements from an array.

TEST(ArrayOperationsTest, SingleElementAccess)
{

    using StringUtilities::dtos;

    arma::fmat array = {
        { 1,  2,  3,  4,  5},
        { 6,  7,  8,  9, 10},
        {11, 12, 13, 14, 15}
    };

    ASSERT_EQ(11, array(2, 0));
    ASSERT_EQ(3, array(0, 2));

    ASSERT_EQ(3, array.n_rows);
    ASSERT_EQ(5, array.n_cols);

}








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
    arma::fmat tinyArray = {
        {0.1, 0.2, 0.3},
        {0.4, 0.5, 0.6},
        {0.7, 0.8, 0.9}
    };

    // Expected result after rotating the smallArray by 90 degrees
    arma::fmat rotatedTinyArray90 = {
        {0.7, 0.4, 0.1},
        {0.8, 0.5, 0.2},
        {0.9, 0.6, 0.3}
    };

    // Test array where there are no zero (0) values at the sides.
    arma::fmat smallArray = {
        {1., 2., 3.},
        {4., 5., 6.},
        {7., 8., 9.}
    };

    // Expected result after rotating the smallArray by 90 degrees
    arma::fmat rotatedSmallArray90 = {
        {7., 4., 1.},
        {8., 5., 2.},
        {9., 6., 3.}
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
        {0. , 7., 4., 1., 0.},
        {0. , 8., 5., 2., 0.},
        {0. , 9., 6., 3., 0.},
        {0. , 0., 0., 0., 0.}
    };


    
    arma::fmat arr = ArrayOperations::rotateArray(tinyArray, deg2rad(90.0));
    checkArraysToBeEqual(rotatedTinyArray90, arr);

//    printArray(smallArray, "Original smallArray");
//    printArray(arr, "Rotated smallArray");

    arr = ArrayOperations::rotateArray(tinyArray, deg2rad(0.0));
    checkArraysToBeEqual(tinyArray, arr);

//    printArray(smallArray, "Original smallArray");
//    printArray(arr, "Rotated smallArray");
    




    arr = ArrayOperations::rotateArray(smallArray, deg2rad(90.0));
    checkArraysToBeEqual(rotatedSmallArray90, arr);

//    printArray(smallArray, "Original smallArray");
//    printArray(arr, "Rotated smallArray");

    arr = ArrayOperations::rotateArray(smallArray, deg2rad(0.0));
    checkArraysToBeEqual(smallArray, arr);

//    printArray(smallArray, "Original smallArray");
//    printArray(arr, "Rotated smallArray");
    




    arr = ArrayOperations::rotateArray(smallNulledArray, deg2rad(90.0));
    checkArraysToBeEqual(rotatedSmallNulledArray90, arr);

//    printArray(smallNulledArray, "Original smallNulledArray");
//    printArray(arr, "Rotated smallNulledArray");




    arr = ArrayOperations::rotateArray(smallNulledArray, deg2rad(0.0));
    checkArraysToBeEqual(smallNulledArray, arr);

//    printArray(smallNulledArray, "Original smallNulledArray");
//    printArray(arr, "Rotated smallNulledArray");




    arr = ArrayOperations::rotateArray(bigNulledArray, deg2rad(90.0));
    EXPECT_FLOAT_EQ(50.0, arr(9, 9));
    EXPECT_FLOAT_EQ(10.0, arr(5, 5));
    EXPECT_FLOAT_EQ(10.0, arr(13, 13));
    EXPECT_FLOAT_EQ(40.0, arr(13, 5));
    EXPECT_FLOAT_EQ(20.0, arr(5, 13));

    arr = ArrayOperations::rotateArray(bigNulledArray, deg2rad(0.0));
    checkArraysToBeEqual(bigNulledArray, arr);


//    printArray(bigNulledArray, "Original bigNulledArray");
//    printArray(arr, "Rotated bigNulledArray");

}


// Test what the effect is of small rebinning actions

TEST(ArrayOperationsTest, DISABLED_RebinningEffect)
{

    arma::fmat array(1024, 1024, arma::fill::eye);

    arma::fmat rebinnedArray = rebin(array, 128, 1, 512, 512);
    printArray(rebinnedArray, "ArrayOperationsTest.RebinningEffect");

    EXPECT_EQ(9, rebinnedArray.n_rows);
    EXPECT_EQ(9, rebinnedArray.n_cols);
    EXPECT_EQ(1024, arma::accu(rebinnedArray));

    rebinnedArray = rebin(array, 8, 8);
    printArray(rebinnedArray, "ArrayOperationsTest.RebinningEffect");

    EXPECT_EQ(8, rebinnedArray.n_rows);
    EXPECT_EQ(8, rebinnedArray.n_cols);
    EXPECT_EQ(1024, arma::accu(rebinnedArray));

}

TEST(ArrayOperationsTest, Rebinning)
{
    arma::fmat array(10, 10, arma::fill::eye);
    arma::fmat testArray(5, 5, arma::fill::eye);
    testArray *= 2;

    arma::fmat rebinnedArray = rebin(array, 5, 5);
    checkArraysToBeEqual(rebinnedArray, testArray);

    array = {
        {1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5},
        {1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5},
        {1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5},

        {2, 2, 2, 1, 1, 1, 4, 4, 4, 5, 5, 5, 4, 4, 4},
        {2, 2, 2, 1, 1, 1, 4, 4, 4, 5, 5, 5, 4, 4, 4},
        {2, 2, 2, 1, 1, 1, 4, 4, 4, 5, 5, 5, 4, 4, 4},

        {3, 3, 3, 4, 4, 4, 1, 1, 1, 4, 4, 4, 3, 3, 3},
        {3, 3, 3, 4, 4, 4, 1, 1, 1, 4, 4, 4, 3, 3, 3},
        {3, 3, 3, 4, 4, 4, 1, 1, 1, 4, 4, 4, 3, 3, 3},

        {4, 4, 4, 5, 5, 5, 4, 4, 4, 1, 1, 1, 2, 2, 2},
        {4, 4, 4, 5, 5, 5, 4, 4, 4, 1, 1, 1, 2, 2, 2},
        {4, 4, 4, 5, 5, 5, 4, 4, 4, 1, 1, 1, 2, 2, 2},

        {5, 5, 5, 4, 4, 4, 3, 3, 3, 2, 2, 2, 1, 1, 1},
        {5, 5, 5, 4, 4, 4, 3, 3, 3, 2, 2, 2, 1, 1, 1},
        {5, 5, 5, 4, 4, 4, 3, 3, 3, 2, 2, 2, 1, 1, 1}
    };

    testArray = {
        {1 * 9, 2 * 9, 3 * 9, 4 * 9, 5 * 9},
        {2 * 9, 1 * 9, 4 * 9, 5 * 9, 4 * 9},
        {3 * 9, 4 * 9, 1 * 9, 4 * 9, 3 * 9},
        {4 * 9, 5 * 9, 4 * 9, 1 * 9, 2 * 9},
        {5 * 9, 4 * 9, 3 * 9, 2 * 9, 1 * 9}
    };

    rebinnedArray = rebin(array, 5, 5);
    checkArraysToBeEqual(rebinnedArray, testArray);


    testArray = {
        { 37,  81, 113},
        { 81,  69,  81},
        {113,  81,  37},
    };

    rebinnedArray = rebin(array, 3, 3);
    checkArraysToBeEqual(rebinnedArray, testArray);


    // Only square arrays are currently supported
    ASSERT_THROW(rebinnedArray = rebin(array, 3, 5), IllegalArgumentException);

    // Only down scaling is currently supported
    ASSERT_THROW(rebinnedArray = rebin(array, 17, 5), IllegalArgumentException);

    // Rebinned array must be integer fraction of original array
    ASSERT_THROW(rebinnedArray = rebin(array, 4, 4), IllegalArgumentException);
}

TEST(ArrayOperationsTest, RebinningCenter)
{
    arma::fmat sourceArray;
    arma::fmat testArray;
    arma::fmat rebinnedArray;


    // This is the simplest case where the array
    //   * is square and of odd size so
    //   * the center element is perfectly in the middle

    sourceArray = {  // 9 x 9, center element is (4, 4)
        {0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 1, 1, 0, 0, 0},
        {0, 0, 0, 1, 1, 1, 0, 0, 0},
        {0, 0, 0, 1, 1, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0},
    };

    rebinnedArray = rebin(sourceArray, 3, 1, 4, 4);

    testArray = {
        {0, 0, 0},
        {0, 9, 0},
        {0, 0, 0},
    };

    checkArraysToBeEqual(rebinnedArray, testArray);

}




TEST(ArrayOperationsTest, RebinningGrow)
{
    arma::fmat sourceArray;
    arma::fmat testArray;
    arma::fmat rebinnedArray;


    sourceArray = {
        {0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 1, 2, 3, 0, 0, 0, 0},
        {0, 0, 4, 5, 6, 0, 0, 0, 2},
        {0, 0, 7, 8, 9, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 3},
        {1, 0, 0, 2, 0, 0, 3, 0, 0}
    };

    testArray = {
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1},
        {0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 4, 4, 5, 5, 6, 6, 0, 0, 0, 0, 0, 0, 2, 2},
        {0, 0, 0, 0, 4, 4, 5, 5, 6, 6, 0, 0, 0, 0, 0, 0, 2, 2},
        {0, 0, 0, 0, 7, 7, 8, 8, 9, 9, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 7, 7, 8, 8, 9, 9, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3},
        {1, 1, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0},
        {1, 1, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0}
    };

    ASSERT_THROW(rebinnedArray = rebin(sourceArray, 1, 2, 2, 3), UnsupportedException);

    // checkArraysToBeEqual(rebinnedArray, testArray);

}






TEST(ArrayOperationsTest, RebinningCenterOdd)
{
    arma::fmat sourceArray;
    arma::fmat testArray;
    arma::fmat rebinnedArray;


    sourceArray = {  // 10 x 10, center element is (4, 4)
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 1, 1, 0, 0, 0, 0},
        {0, 0, 0, 1, 1, 1, 0, 0, 0, 2},
        {0, 0, 0, 1, 1, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 3},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 2, 0, 0, 3, 0, 0},
    };

    testArray = {
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1},
        {0, 0, 9, 0, 2},
        {0, 0, 0, 0, 3},
        {0, 1, 2, 3, 0},
    };

    rebinnedArray = rebin(sourceArray, 3, 1, 4, 4);
    checkArraysToBeEqual(rebinnedArray, testArray);

    testArray = {
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 3, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    };

    rebinnedArray = rebin(sourceArray, 3, 3, 4, 4);
    checkArraysToBeEqual(rebinnedArray, testArray);

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