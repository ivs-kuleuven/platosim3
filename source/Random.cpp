#include "Random.h"


// Constructor
// INPUT:
//    location: real
//    scale: real > 0
//    shape: real

skew_normal_distribution::skew_normal_distribution(double location, double scale, double shape)
: location(location), scale(scale), shape(shape)
{
    assert(scale > 0.0);
    standardNormal = normal_distribution<double>(0.0, 1.0);
}




// Drawing random number


double skew_normal_distribution::operator()(mt19937 &engine)
{
   const double u = standardNormal(engine);
   const double v = standardNormal(engine);

   double randomNumber = shape/sqrt(1+shape*shape) * abs(u) + 1/sqrt(1+shape*shape) * v;
   randomNumber = randomNumber * scale + location;
   return randomNumber;
}

