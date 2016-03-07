#ifndef POLYNOMIAL_1D_H
#define POLYNOMIAL_1D_H


#include <iostream>

#include "Exceptions.h"
#include "Logger.h"
#include "StringUtilities.h"

class Polynomial1D 
{

    public:
        Polynomial1D();
        Polynomial1D(int degree, vector<double> coefficients);
        ~Polynomial1D();

        double getCoefficient(int index);
        double operator()(double x);    

    private:
        int degree;
        vector<double> coefficients;
};


#endif /* POLYNOMIAL_1D_H */
