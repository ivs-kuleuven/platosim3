#include <iostream>
#include <iomanip>


#include "Units.h"
#include "Logger.h"

#include "ArrayOperations.h"


using std::cout;
using std::endl;
using std::stringstream;
using std::setprecision;
using std::setiosflags;

// Functions used locally

arma::Mat<float> rotateArrayNearestNeighbor(arma::Mat<float> arr, double angle);
arma::Mat<float> rotateArrayBilinear(arma::Mat<float> arr, double angle);







/**
 * \brief      Rotate a 2D Armadillo array over an arbitrairy angle
 * 
 * \details
 * 
 * This function currently uses a simple (not optimized) version of bilinear interpolation
 * to find the proper pixel value.
 * 
 * The rotation is done with respect to the positive x-axis and counter clockwise.
 *
 * \param[in]  arr    a 2D array of floats
 * \param[in]  angle  the angle over which the array must be rotated [degrees]
 *
 * \return     the rotated array with the same dimensions as the original array
 */
arma::Mat<float> rotateArray(arma::Mat<float> arr, double angle)
{
    return rotateArrayBilinear(arr, angle);
}







/**
 * \brief      Print the 2D array to std::cout
 *
 * \param[in]  arr   a 2D armadillo array
 * \param[in]  msg   a message that is printed before the array is printed
 */
void printArray(arma::Mat<float> arr, string msg)
{
    stringstream fullMessage;

    fullMessage << setiosflags(ios::fixed);
    fullMessage << setprecision(4);

    fullMessage << msg << endl;
    
    for (int row = 0; row < arr.n_rows; row++)
    {
        for (int column = 0; column < arr.n_cols; column++)
        {
            fullMessage << arr(row, column) << " ";
        }
        fullMessage << endl;
    }
    fullMessage << "Sum = " << arma::accu(arr) << endl;

    cout << fullMessage.str();
}






/**
 * \brief      Rotate a 2D armadillo array and use nearest neighbor selection 
 *
 * \param[in]  origArray  a 2D armadillo array of floats
 * \param[in]  angle      the degrees by which to rotate the array [degrees]
 *
 * \return     a rotated 2D armadillo array with the same dimensions
 */
arma::Mat<float> rotateArrayNearestNeighbor(arma::Mat<float> origArray, double angle)
{
    arma::Mat<float> rotatedArray(arma::size(origArray), arma::fill::zeros);

    int width = origArray.n_rows;
    int height = origArray.n_cols;

    float background = 0.0;  // this is the background used when no source pixel is matching
    float rads = deg2rad(angle);
    float cs = std::cos(-rads);   // precalculate these values, we negate the angle because we are going back from destination to source
    float ss = std::sin(-rads);
    float xcenter = ((float)(width) / 2.0);
    float ycenter = ((float)(height) / 2.0);

    Log.debug("xcenter, ycenter = " + to_string(xcenter) + ", " + to_string(ycenter));

    for (int row = 0; row < width; ++row)
    {
       for (int column = 0; column < height; ++column)
       {
          // Find the pixel of the original image that is rotated to (row, column)
          // The rotation formula assumes that origin = top-left and y points down
          float xPrime = ycenter + ((float)(row)+0.5-ycenter)*cs - ((float)(column)+0.5-xcenter)*ss;
          float yPrime = xcenter + ((float)(row)+0.5-ycenter)*ss + ((float)(column)+0.5-xcenter)*cs;

          // Find the pixel coordinates of the source pixel
          int xPixel = std::round(xPrime-0.5);
          int yPixel = std::round(yPrime-0.5);

          // Log.debug("row, col | xPrime, yPrime | xPixel, yPixel = " + to_string(row) + ", " + to_string(column) + " | " + to_string(xPrime) + ", " + to_string(yPrime) + " | " + to_string(xPixel) + ", " + to_string(yPixel));

          float pixel = background; // in case there is no original pixel
          
          // Get the pixel value from the source if the pixel falls within the original array

          if (xPixel >= 0 && xPixel < width && yPixel >= 0 && yPixel < height)
          {
             pixel = origArray(xPixel, yPixel);
          }
          rotatedArray(row, column) = pixel;
       }
    }

    return rotatedArray;
}








/**
 * \brief      Rotate a 2D armadillo array and use bilinear interpolation to determine the pixel value 
 *
 * \param[in]  origArray  a 2D armadillo array of floats
 * \param[in]  angle      the degrees by which to rotate the array [degrees]
 *
 * \return     a rotated 2D armadillo array with the same dimensions
 */
arma::Mat<float> rotateArrayBilinear(arma::Mat<float> origArray, double angle)
{
    arma::Mat<float> rotatedArray(arma::size(origArray), arma::fill::zeros);

    int width = origArray.n_rows;
    int height = origArray.n_cols;
    
    double rads = deg2rad(angle);
    double cs = std::cos(-rads);  // precalculate these values
    double ss = std::sin(-rads);

    double cX = (double)width / 2.0;
    double cY = (double)height / 2.0;

    int q11x, q12x, q21x, q22x;
    int q11y, q12y, q21y, q22y;

    Log.debug("cX, cY = " + to_string(cX) + ", " + to_string(cY));

    for (int row=0; row < width; row++)
    {
        double relX = (double)row+0.5 - cY;
        for (int column=0; column < height; column++)
        {
            double relY = (double)column+0.5 - cX;

            double xPrime = cY + relX*cs - relY*ss - 0.5;
            double yPrime = cX + relX*ss + relY*cs - 0.5;

            // Log.debug("row, col | xPrime, yPrime = " + to_string(row) + ", " + to_string(column) + " | " + to_string(xPrime) + ", " + to_string(yPrime));

            int xPixel = std::round(xPrime);
            int yPixel = std::round(yPrime);

            // Log.debug("xPixel, yPixel = " + to_string(xPixel) + ", " + to_string(yPixel));

            if (yPrime <= (double)yPixel)
            {
                //Log.debug("yPrime in lower part of yPixel");
                q12y = q22y = yPixel - 1;
                q11y = q21y = yPixel;
            }
            else 
            {
                q12y = q22y = yPixel;
                q11y = q21y = yPixel + 1;
            }

            if (xPrime <= (double)xPixel)
            {
                //Log.debug("xPrime in lower part of xPixel");
                q12x = q11x = xPixel - 1;
                q22x = q21x = xPixel;
            }
            else
            {
                q12x = q11x = xPixel;
                q22x = q21x = xPixel + 1;
            }

            // Log.debug("Q12x, Q12y = " + to_string(q12x) + ", " + to_string(q12y));
            // Log.debug("Q22x, Q22y = " + to_string(q22x) + ", " + to_string(q22y));
            // Log.debug("Q11x, Q11y = " + to_string(q11x) + ", " + to_string(q11y));
            // Log.debug("Q21x, Q21y = " + to_string(q21x) + ", " + to_string(q21y));

            double factor1;
            double factor2;

            unsigned int q11, q12, q21, q22;

            // We need to get the four nearest neighbooring pixels.
            // Pixels which are past the border of the image are clamped to the border already.
            if (q11x < 0 or q11x > width-1 or q11y < 0 or q11y > height-1)
                q11 = 0;
            else
                q11 = origArray(q11x, q11y);
            if (q12x < 0 or q12x > width-1 or q12y < 0 or q12y > height-1)
                q12 = 0;
            else
                q12 = origArray(q12x, q12y);
            if (q21x < 0 or q21x > width-1 or q21y < 0 or q21y > height-1)
                q21 = 0;
            else
                q21 = origArray(q21x, q21y);
            if (q22x < 0 or q22x > width-1 or q22y < 0 or q22y > height-1)
                q22 = 0;
            else
                q22 = origArray(q22x, q22y);      

            if ( q21x == q11x ) // special case to avoid division by zero
            {
                factor1 = 1; // the same X coordinate, so just force the calculation to one point
                factor2 = 0;
            }
            else
            {
                factor1 = (((double)q21x - xPrime)/((double)q21x - (double)q11x));
                factor2 = ((xPrime - (double)q11x)/((double)q21x - (double)q11x));
            }
            double R1 = factor1 * (double)q11 + factor2*(double)q21;
            double R2 = factor1 * (double)q12 + factor2*(double)q22;

            double factor3;
            double factor4;

            if (q12y == q11y) // special case to avoid division by zero
            {
                factor3 = 1;
                factor4 = 0;
            }
            else
            {
                factor3 = ((double) q12y - yPrime)/((double)q12y - (double)q11y);
                factor4 = (yPrime - (double)q11y)/((double)q12y - (double)q11y);
            }
            
            double finalR = (factor3 * R1) + (factor4 * R2);

            rotatedArray(row, column) = finalR;
        }
    }

    return rotatedArray;
}





