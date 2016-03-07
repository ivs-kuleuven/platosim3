#ifndef ARRAY_OPERATIONS_H
#define ARRAY_OPERATIONS_H

#include <string>
#include <iostream>
#include <iomanip>

#include "armadillo"

#include "Units.h"
#include "Logger.h"


using namespace std;


arma::Mat<float> rotateArray(arma::Mat<float> arr, double angle);
void printArray(arma::Mat<float> arr, string msg);



#endif /* ARRAY_OPERATIONS_H */

