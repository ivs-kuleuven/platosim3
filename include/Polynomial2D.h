#ifndef POLYNOMIAL_2D_H
#define POLYNOMIAL_2D_H


#include <iostream>

#include "Exceptions.h"
#include "Logger.h"
#include "StringUtilities.h"

class Polynomial2D 
{

    public:
        Polynomial2D();
        Polynomial2D(int degree, vector<double> coefficients);
        ~Polynomial2D();

        double getCoefficient(int index);
        double operator()(double x, double y);    

    private:
        int degree;
        vector<double> coefficients;
};


#endif /* POLYNOMIAL_2D_H */
