#include "Polynomial1D.h"

using namespace std;






/**
 * \brief      Constructor for the Polynomial1D class
 * 
 * \detail
 * 
 * The general 1D polynomial of degree n is
 * 
 * \f[P(x) = c_{0} + c_{1} x + c_{2} x^{2} + ... + c_{n} x^{n}\f]
 *
 * \param[in]  deg    the degree of the Polynomial
 * \param      coeff  the coeficients of the polynomial
 * 
 */
Polynomial1D::Polynomial1D(int deg, vector<double> coeff)
{
    degree = deg;

    // Copying the argument coefficients
    coefficients = coeff;
}







/**
 * \brief      Destructor, free coefficients memory
 */
Polynomial1D::~Polynomial1D()
{
}








/**
 * \brief      Evaluate the 1D polynomial in x
 *
 * \param[in]  x     x-value
 *
 * \return     y the evaluated polynomial
 */
double Polynomial1D::operator()(double x)
{

    using StringUtilities::dtos;

    double sum = coefficients[0];
    double xPow = 1.0;

    for(unsigned int i = 1; i <= degree; i++)
    {
        xPow *= x;
        sum += xPow * coefficients[i];
    }

    return sum;
}

