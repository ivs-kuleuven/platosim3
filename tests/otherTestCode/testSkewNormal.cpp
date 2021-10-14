
// Code that can be used to test the skew Normal stuff
//
// Compile with
// clang++ -o testSkewNormal testSkewNormal.cpp ../../source/Random.cpp  -I ../../include -stdlib=libc++ -std=c++14
//

#include <iostream>
#include <random>

#include "Random.h"


using namespace std;


int main()
{
    mt19937 generator;
    generator.seed(12345678);

    const double location = 2500.0;
    const double scale = 2000.0;
    const double shape = 30.0;

    skew_normal_distribution distrib(location, scale, shape);

    for (unsigned int n = 0; n < 100000; n++)
    {
        cout << distrib(generator) << endl;
    }

    return 0;
}
