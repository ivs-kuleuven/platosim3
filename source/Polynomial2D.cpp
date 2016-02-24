#include "Polynomial2D.h"

using namespace std;



/**
 * \brief      Calculate the number of terms for the 2D polynomial of degree n
 *
 * \param[in]  n     the degree of the polynomial (n <= 3)
 *
 * \return     the number of terms
 */
static unsigned int nTerms(unsigned int n)
{
    unsigned int ret = 0;
    for(unsigned int i = 1; i <= n+1; ++i)
        ret += i;
    return ret;
}




Polynomial2D::Polynomial2D(int deg, double coeff[])
{
    degree = deg;
    unsigned int n_terms = nTerms(degree);

    //Log.debug("Number of terms for " + to_string(degree) + " is " + to_string(n_terms));

    coefficients = new double[n_terms];

    if (degree > 3)
    {
        throw UnsupportedException("Polynomial2D::Constructor: Evaluation of a 2D polynomial of degree > 3 is not yet supported.");
    }

    for( unsigned int i = 0; i < n_terms; i++)
    {
        coefficients[i] = coeff[i];
    }
}

Polynomial2D::~Polynomial2D()
{
    delete [] coefficients;
}

double Polynomial2D::evaluateAt(double x, double y)
{

    using StringUtilities::dtos;

    double sum = coefficients[0];
    double xPow = 1.0;
    double yPow = 1.0;

    Log.debug("sum = " + dtos(sum));

    for(unsigned int i = 1; i <= degree; i++)
    {
        xPow *= x;
        yPow *= y;
        sum += xPow * coefficients[i] + yPow * coefficients[degree+i];
    }

    // now we take care of the mixing terms
    // this code can probably be optimized using Pascal's triangle
    // but we now do this manually

    // no mixing terms for degree == 1

    unsigned int idx = degree * 2 + 1;
    
    if (degree >= 2)
    {
        // mixing term for degrees == 2 is c11 * x * y
        sum += coefficients[idx++] * x * y;
    }

    if (degree >= 3)
    {
        // mixing terms for degree == 3 are c12 * x * y^2 and c21 * x^2 * y
        sum += coefficients[idx++] * x * y * y;
        sum += coefficients[idx++] * x * x * y;
    }

    return sum;
}

