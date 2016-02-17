#include "gtest/gtest.h"

#include "Units.h"
#include "Constants.h"
#include "ArrayOperations.h"

using namespace std;


void checkArraysToBeEqual(arma::Mat<float> arr1, arma::Mat<float> arr2);



TEST(ArrayOperationsTest, Rotation)
{

    arma::Mat<float> bigNulledArray = {
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

    arma::Mat<float> smallArray = {
        {1., 2., 3.},
        {4., 5., 6.},
        {7., 8., 9.}
    };

    arma::Mat<float> rotatedSmallArray90 = {
        {3., 6., 9.},
        {2., 5., 8.},
        {1., 4., 7.}
    };

    arma::Mat<float> smallNulledArray = {
        {0. , 0., 0., 0., 0.},
        {0. , 1., 2., 3., 0.},
        {0. , 4., 5., 6., 0.},
        {0. , 7., 8., 9., 0.},
        {0. , 0., 0., 0., 0.}
    };

    arma::Mat<float> rotatedSmallNulledArray90 = {
        {0. , 0., 0., 0., 0.},
        {0. , 3., 6., 9., 0.},
        {0. , 2., 5., 8., 0.},
        {0. , 1., 4., 7., 0.},
        {0. , 0., 0., 0., 0.}
    };


    
    arma::Mat<float> arr = rotateArray(smallArray, 90.0);
    checkArraysToBeEqual(rotatedSmallArray90, arr);

    //printArray(smallArray, "Original smallArray");
    //printArray(arr, "Rotated smallArray");
    




    arr = rotateArray(smallNulledArray, 90.0);
    checkArraysToBeEqual(rotatedSmallNulledArray90, arr);

    //printArray(smallNulledArray, "Original smallNulledArray");
    //printArray(arr, "Rotated smallNulledArray");




    arr = rotateArray(bigNulledArray, 90.0);

    //printArray(bigNulledArray, "Original bigNulledArray");
    //printArray(arr, "Rotated bigNulledArray");

}



void checkArraysToBeEqual(arma::Mat<float> arr1, arma::Mat<float> arr2)
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