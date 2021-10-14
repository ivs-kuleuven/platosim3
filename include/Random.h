#ifndef RANDOM_H
#define RANDOM_H

#include <cassert>
#include <cmath>
#include <iostream>
#include <random>

using namespace std;

// A random variable Z has a skew-normal distribution with shape parameter lambda, denoted by Z ~ SN(lambda), 
// if its density is given by f(z, lambda) = 2 * Phi(lambda*z) * phi(z) where Phi() and phi() are the
// standard normal cumulative distribution function and the standard normal probability. 

class skew_normal_distribution 
{
    public:

        skew_normal_distribution(double location=0.0, double scale=1.0, double shape=0.0);
        double operator()(mt19937 &engine);

    protected:

    private:

        double location;
        double scale;
        double shape;

        normal_distribution<double> standardNormal;        

}; // end class skew_normal_distribution


#endif
