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




/**
 * @brief      Default constructor, x, y = 1.0
 */
Polynomial2D::Polynomial2D()
: degree(0), coefficients({1.0})
{

}




/**
 * \brief      Constructor for the Polynomial2D class
 * 
 * \detail
 * 
 * The general 2D polynomial of degree n is
 * 
 * \f[P(x, y) = c_{00} + c_{10} x + ... + c_{n0} x^{n} + c_{01} y 
 *     + ... + c_{0n} y^{n} + c_{11} x y + c_{12} x y^{2} 
 *     + ... + c_{1(n−1)} x y^{(n−1)} + ... + c_{(n−1)1} x^{(n−1)} y
 * \f]
 *
 * \param[in]  deg    the degree of the Polynomial (n <= 3)
 * \param      coeff  the coeficients of the polynomial
 * 
 * \exception UnsupportedException for degrees higher than 3
 */
Polynomial2D::Polynomial2D(int deg, vector<double> coeff)
{
    degree = deg;

    unsigned int n_terms = nTerms(degree);

    if (degree > 3)
    {
        throw UnsupportedException("Polynomial2D::Constructor: Evaluation of a 2D polynomial of degree > 3 is not yet supported.");
    }

    if (coeff.size() != n_terms)
    {
        throw IllegalArgumentException("Polynomial2D::Constructor: Expected " + to_string(n_terms) + " coefficients for degree " + to_string(degree) + ", got " + to_string(coeff.size()) + " coefficients.");
    }

    // Copying the argument coefficients
    coefficients = coeff;
}







/**
 * \brief      Destructor, free coefficients memory
 */
Polynomial2D::~Polynomial2D()
{
}








/**
 * \brief      Evaluate the 2D polynomial in x and y
 *
 * \param[in]  x     x-value
 * \param[in]  y     y-value
 *
 * \return     z the evaluated polynomial
 */
double Polynomial2D::operator()(double x, double y)
{

    using StringUtilities::dtos;

    double sum = coefficients[0];
    double xPow = 1.0;
    double yPow = 1.0;

    //Log.debug("sum = " + dtos(sum));

    for(unsigned int i = 1; i <= degree; i++)
    {
        xPow *= x;
        yPow *= y;
        sum += xPow * coefficients[i] + yPow * coefficients[degree+i];
        //Log.debug("sum = " + dtos(sum));
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
        //Log.debug("sum = " + dtos(sum));
    }

    if (degree >= 3)
    {
        // mixing terms for degree == 3 are c12 * x * y^2 and c21 * x^2 * y
        sum += coefficients[idx++] * x * y * y;
        //Log.debug("sum = " + dtos(sum));

        sum += coefficients[idx++] * x * x * y;
        //Log.debug("sum = " + dtos(sum));
    }

    return sum;
}

