#include "ArrayOperations.h"


using namespace std;




namespace ArrayOperations
{

// Functions Prototypes 

arma::fmat rotateArrayNearestNeighbor(arma::fmat arr, double angle);
arma::fmat rotateArrayBilinear(arma::fmat arr, double angle);




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
 * \param[in]  angle  the angle over which the array must be rotated [radians]
 *
 * \return     the rotated array with the same dimensions as the original array
 */
arma::fmat rotateArray(arma::fmat arr, double angle)
{
    return rotateArrayBilinear(arr, angle);
}







/**
 * \brief      Print the 2D array to std::cout
 *
 * \param[in]  arr   a 2D armadillo array
 * \param[in]  msg   a message that is printed before the array is printed
 */
void printArray(arma::fmat arr, string msg)
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
 * \param[in]  angle      the degrees by which to rotate the array [radians]
 *
 * \return     a rotated 2D armadillo array with the same dimensions
 */
arma::fmat rotateArrayNearestNeighbor(arma::fmat origArray, double angle)
{
    arma::fmat rotatedArray(arma::size(origArray), arma::fill::zeros);

    int width = origArray.n_rows;
    int height = origArray.n_cols;

    float background = 0.0;  // this is the background used when no source pixel is matching
    
    // precalculate these values, 
    // we should negate the angle because we are going back from destination to source, but
    // since the rotation is on the positive x-axis going counter-clockwise, and the array
    // is flipped vertically when displayed,  we have to negate the angle again which results
    // is just the angle.
    float cs = std::cos(angle);
    float ss = std::sin(angle);
    float xcenter = ((float)(width) / 2.0);
    float ycenter = ((float)(height) / 2.0);

    // Log.debug("xcenter, ycenter = " + to_string(xcenter) + ", " + to_string(ycenter));

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
 * \param[in]  angle      the degrees by which to rotate the array [radians]
 *
 * \return     a rotated 2D armadillo array with the same dimensions
 */
arma::fmat rotateArrayBilinear(arma::fmat origArray, double angle)
{
    arma::fmat rotatedArray(arma::size(origArray), arma::fill::zeros);

    int width = origArray.n_rows;
    int height = origArray.n_cols;
    
    double cs = std::cos(angle);  // precalculate these values
    double ss = std::sin(angle);

    // Log.debug("rotateArrayBilinear: cs, ss = " + to_string(cs) + ", " + to_string(ss));

    double cX = (double)width / 2.0;
    double cY = (double)height / 2.0;

    int q11x, q12x, q21x, q22x;
    int q11y, q12y, q21y, q22y;

    // Log.debug("rotateArrayBilinear: cX, cY = " + to_string(cX) + ", " + to_string(cY));

    for (int row=0; row < width; row++)
    {
        double relX = (double)row+0.5 - cY;
        for (int column=0; column < height; column++)
        {
            double relY = (double)column+0.5 - cX;

            double xPrime = cY + relX*cs - relY*ss - 0.5;
            double yPrime = cX + relX*ss + relY*cs - 0.5;

            // Log.debug("rotateArrayBilinear: row, col | xPrime, yPrime = " + to_string(row) + ", " + to_string(column) + " | " + to_string(xPrime) + ", " + to_string(yPrime));

            int xPixel = std::round(xPrime);
            int yPixel = std::round(yPrime);

            // Log.debug("rotateArrayBilinear: xPixel, yPixel = " + to_string(xPixel) + ", " + to_string(yPixel));

            if (yPrime <= (double)yPixel)
            {
                //Log.debug("rotateArrayBilinear: yPrime in lower part of yPixel");
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
                //Log.debug("rotateArrayBilinear: xPrime in lower part of xPixel");
                q12x = q11x = xPixel - 1;
                q22x = q21x = xPixel;
            }
            else
            {
                q12x = q11x = xPixel;
                q22x = q21x = xPixel + 1;
            }

            // Log.debug("rotateArrayBilinear: Q12x, Q12y = " + to_string(q12x) + ", " + to_string(q12y));
            // Log.debug("rotateArrayBilinear: Q22x, Q22y = " + to_string(q22x) + ", " + to_string(q22y));
            // Log.debug("rotateArrayBilinear: Q11x, Q11y = " + to_string(q11x) + ", " + to_string(q11y));
            // Log.debug("rotateArrayBilinear: Q21x, Q21y = " + to_string(q21x) + ", " + to_string(q21y));

            double factor1;
            double factor2;

            double q11, q12, q21, q22;

            // We need to get the four nearest neighbooring pixels.
            // Pixels which are past the border of the image are clamped to the border already.
            if (q11x < 0 or q11x > width-1 or q11y < 0 or q11y > height-1)
                q11 = 0.0;
            else
                q11 = origArray(q11x, q11y);
            if (q12x < 0 or q12x > width-1 or q12y < 0 or q12y > height-1)
                q12 = 0.0;
            else
                q12 = origArray(q12x, q12y);
            if (q21x < 0 or q21x > width-1 or q21y < 0 or q21y > height-1)
                q21 = 0.0;
            else
                q21 = origArray(q21x, q21y);
            if (q22x < 0 or q22x > width-1 or q22y < 0 or q22y > height-1)
                q22 = 0.0;
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

            // Log.debug("rotateArrayBilinear: factor1, factor2 = " + to_string(factor1) + ", " + to_string(factor2));
            
            double R1 = factor1 * q11 + factor2 * q21;
            double R2 = factor1 * q12 + factor2 * q22;

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

            // Log.debug("rotateArrayBilinear: factor3, factor4 = " + to_string(factor3) + ", " + to_string(factor3));
            
            double finalR = (factor3 * R1) + (factor4 * R2);

            // Log.debug("rotateArrayBilinear: finalR = " + to_string(finalR));

            rotatedArray(row, column) = finalR;
        }
    }

    return rotatedArray;
}








/**
 * @brief      Method that rebins an array to the given new resolution.
 *
 * @details    
 * 
 * The array is considered square and the new dimension must be a multiple or 
 * a factor of the dimension of array.
 * 
 * 
 * 
 * Remember the arrays are column-major.
 *
 * @param[in]  array             the original array that needs rebinning
 * @param[in]  sourceResolution  the number of elements that make up a group in the source
 * @param[in]  targetResolution  the number of elements that make up a group in the target
 * @param[in]  xCenter           the center x point of the original array (zero-based)
 * @param[in]  yCenter           the center y point of the original array (zero-based)
 *
 * @return     The new rebinned array
 */
arma::fmat rebin(arma::fmat array, unsigned int sourceResolution, unsigned int targetResolution,
                       unsigned int xCenter, unsigned int yCenter)
{
    // Create an intermediate array of odd size for which the centre point of the
    // original array lies in the middle of the intermediate array.
    // The new intermediate array shall be big enough so the complete original array fits in.

    unsigned int  sizeIntermediateArray = max(
        (unsigned int) max(xCenter + 1, yCenter + 1),
        (unsigned int) max(array.n_rows - xCenter, array.n_cols - yCenter)) * 2 - 1;

    // The dimension of the intermediate array shall be a multiple of the source resolution.

    sizeIntermediateArray = sizeIntermediateArray + (sizeIntermediateArray % sourceResolution) * 2;

    arma::fmat intermediateArray(sizeIntermediateArray, sizeIntermediateArray, arma::fill::zeros);

    // Copy the original array to the intermediate array and make sure the center point for the original
    // array is now the real center of the intermediate array.

    double centerIntermediateArray = (sizeIntermediateArray - 1) / 2.0;

    intermediateArray(
        arma::span(centerIntermediateArray - xCenter, centerIntermediateArray - xCenter + array.n_rows - 1),
        arma::span(centerIntermediateArray - yCenter, centerIntermediateArray - yCenter + array.n_cols - 1)) =
        array;


    double binning = double(sourceResolution) / double(targetResolution);

    // Size of the rebinned, target array

    unsigned int sizeTargetArray = int(sizeIntermediateArray / binning);

    arma::fmat target;

    if (binning >= 1.0)
    {

        if (sizeTargetArray % 2 == 0)
            sizeTargetArray++;

        target = arma::zeros<arma::fmat>(sizeTargetArray, sizeTargetArray);

        for (int row = 0; row < sizeTargetArray; row++) 
        {
            for (int column = 0; column < sizeTargetArray; column++) 
            {
                if ((row + 1) * binning - 1 < sizeIntermediateArray
                        && (column + 1) * binning - 1 < sizeIntermediateArray)
                {
                    target(row, column) = arma::accu(
                        intermediateArray(arma::span(row * binning, (row + 1) * binning - 1),
                                          arma::span(column * binning, (column + 1) * binning - 1)));
                }
            }
        }
    }
    else
    {
        stringstream fullMessage;

        fullMessage << "ArrayOperations.rebin: Rebinning is only supported to smaller dimensions, i.e. targetResolution (" 
                    << to_string(targetResolution) << "), ";
        fullMessage << "must be smaller than sourceResolution (" << to_string(sourceResolution) << ")." << endl;

        throw UnsupportedException(fullMessage.str());
    }

    return target;
}









/**
 * @brief      Resize an array to the new dimensions.
 *
 * @details
 * 
 * There are currently some restrictions on the rebinning:
 * 
 * \li only square arrays can be rebinned
 * \li rebinning works only towards smaller dimensions
 * \li the dimension of the rebinned array must be an integer fraction of the original array
 * 
 * @param[in]  array       the array that needs to be rebinned
 * @param[in]  n_rows_new  the number of rows in the rebinned array
 * @param[in]  n_cols_new  the number of columns in the new array
 *
 * @return     the rebinned array
 */
arma::fmat rebin(arma::fmat array, unsigned int n_rows_new, unsigned int n_cols_new)
{
    arma::fmat rebinnedArray = arma::zeros<arma::fmat>(n_rows_new, n_cols_new);

    if (n_rows_new != n_cols_new)
    {
        throw IllegalArgumentException(string("ArrayOperations.rebin: ") + 
            "rebinning currently only implemented for a square array, n_rows_new should equal n_cols_new.");
    }

    if (n_rows_new > array.n_rows)
    {
        throw IllegalArgumentException(string("ArrayOperations.rebin: ") + 
            "rebinning currently only implemented to smaller dimensions.");
    }

    if (array.n_rows % n_rows_new)
    {
        throw IllegalArgumentException(string("ArrayOperations.rebin: ") + 
            "new dimensions must be a integer fraction of original dimensions.");

    }

    unsigned int binning = array.n_rows / n_rows_new;

    // Rebinning is simply done by adding all values of the sub-pixels per pixel.

    for (unsigned int row = 0; row < n_rows_new; row++)
    {
        for (unsigned int column = 0; column < n_cols_new; column++)
        {
            const unsigned int beginRow = row * binning;
            const unsigned int beginCol = column * binning;
            const unsigned int endRow = (row + 1) * binning - 1;
            const unsigned int endCol = (column + 1) * binning - 1;

            rebinnedArray(row, column) = arma::accu(array.submat(beginRow, beginCol, endRow, endCol));
        }
    }

    return rebinnedArray;
}

}  // End of namespace ArrayOperations

