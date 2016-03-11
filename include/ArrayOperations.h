#ifndef ARRAY_OPERATIONS_H
#define ARRAY_OPERATIONS_H

#include <string>
#include <iostream>
#include <iomanip>

#include "armadillo"

#include "Units.h"
#include "Logger.h"
#include "Exceptions.h"


using namespace std;


namespace ArrayOperations
{
    arma::Mat<float> rotateArray(arma::Mat<float> arr, double angle);
    void printArray(arma::Mat<float> arr, string msg);
    arma::fmat rebin(arma::fmat array, unsigned int n_rows_new, unsigned int n_cols_new);
    arma::fmat rebin(arma::fmat array, unsigned int sourceResolution, unsigned int targetResolution,
                     unsigned int xCenter, unsigned int yCenter);
}


#endif /* ARRAY_OPERATIONS_H */

