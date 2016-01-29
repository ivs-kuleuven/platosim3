
#ifndef TABULATEDFUNCTION_H
#define TABULATEDFUNCTION_H

#include <vector>
#include <cmath>
#include "Logger.h"

using namespace std;





enum Extrap_Method {Linear_Extrapolation};
enum Interp_Method {Linear_Interpolation, Spline_Interpolation};



template <class Type>
class TabulatedFunction
{
    public:

        TabulatedFunction();
        ~TabulatedFunction();

        void init(Type &a, Type &b, unsigned long N);
        void setExtrapolationMethod(Extrap_Method method);
        void setInterpolationMethod(Interp_Method method);

        double operator()(const double xvalue);

        double integrate(double lower, double upper);
        void setAccuracy(const double accur);

    private:

        Extrap_Method extrap_method;
        Interp_Method interp_method;

        double accuracy;                // fraction accuracy for convergence

        unsigned long Nvalues;          // Number of tabulated values
        double *x;                      // The abscissa values
        double *y;                      // The ordinate values


        double *DDy;                    // Second derivatives (for spline interpolation)
        double *temp;                   // temporary storage (for spline interpolation)


        void init_spline();
        double spline_evaluate(const double xvalue);
        double rational_evaluate(const double xvalue);
        double linear_evaluate(const double xvalue);
};






/**
 *  \brief Default constructor
 *  
 *  \note Initialise the array pointers to 0, so that we can harmlessly delete [] them.
 */

template <class Type>
TabulatedFunction<Type>::TabulatedFunction()
: extrap_method(Linear_Extrapolation), interp_method(Spline_Interpolation),
  accuracy(1.0e-5), Nvalues(0), x(0), y(0), DDy(0), temp(0)
{

} 









/**
 * \brief Destructor
 */

template <class Type>
TabulatedFunction<Type>::~TabulatedFunction()
{
    delete [] x;
    delete [] y;
    delete [] DDy;
    delete [] temp;
}











/**
 * \brief Copy the user given arrays into internal arrays, and
 *        initialize the interpolation procedure.
 *
 * \param a: array with abscissa values of the tabulated function
 * \param b: corresponding array with ordinate values of the tabulated function
 * \param N: only elements 0..N-1 of a[] and b[] will be copied to the internal arrays.
 * 
 * \note: . The abscissa array a must be sorted (ascending), and must be 
 *          free from duplicate values.
 *        . I use a template so that arrays, vectors and valarrays can
 *          be used.
 */

template <class Type>
void TabulatedFunction<Type>::init(Type &a, Type &b, unsigned long N)
{
    // Check whether the given size of the arrays is larger then zero
    // If so, allocate the necessary memory, if not complain.

    if (N > 0)
    {
        // First de-allocate memory. Necessary, if the user wants to
        // use Init multiple times

        Nvalues = N;
        delete [] x; 
        delete [] y;
        delete [] DDy;
        delete [] temp;

        // Then allocate the necessary memory.

        x = new double[N];
        y = new double[N];
        DDy = new double[N];
        temp = new double[N];
    }
    else
    {
        Log.error("TabulatedFunction::init(): size <= 0");
        exit(1);
    }

    // Copy the user given arrays into internal arrays

    for (unsigned long i = 0; i < N; i++)
    {
        x[i] = double(a[i]);
        y[i] = double(b[i]);
    }


    // Check whether the abscissa array is sorted in an ascending way

    for (unsigned long i = 1; i < N; i++)
    {
        if (x[i] <= x[i-1])
        {
            Log.error("TabulatedFunction::init(): x-values not sorted");
            exit(1);
        }
    }

    // Initialize the spline

    init_spline();

}






/**
 * \brief Decide how you wish to compute values out of the x-range, by extrapolation.
 */

template <class Type>
void TabulatedFunction<Type>::setExtrapolationMethod(Extrap_Method method)
{
    if (method == Linear_Extrapolation) 
    {
        extrap_method = method;
    }
    else
    {
        Log.error("TabulatedFunction::setExtrapolationMethod(): Illegal argument");
        exit(1);
    }
}













/**
 * \brief Decide how you wish to compute values inside the x-range, by interpolation.
 */

 template <class Type>
void TabulatedFunction<Type>::setInterpolationMethod(Interp_Method method)
{
    if ((method == Linear_Interpolation) || (method == Spline_Interpolation))
    {
      interp_method = method;
    }
    else
    {
        Log.error("TabulatedFunction::setInterpolationMethod(): Illegal argument");
        exit(1);
    }
}















/**
 * \brief: Evaluate the function by automatic interpolation in the
 *         user given tabulated values.
 */

template <class Type>
double TabulatedFunction<Type>::operator()(const double xvalue)
{
    //  First check whether the class has already been initialized.

    if (Nvalues == 0)
    {
        Log.error("TabulatedFunction::(): TabulatedFunction not initialised.");
        exit(1);
    }

    // Check whether the given value is in its appropriate boundaries
    // and inter/extra-polate accordingly

    if ((xvalue >= x[0]) && (xvalue <= x[Nvalues-1]))
    {
        // Interpolate 

        if (interp_method == Linear_Interpolation)
        {
            return(linear_evaluate(xvalue));
        }

        if (interp_method == Spline_Interpolation)
        {
          return(spline_evaluate(xvalue));
        }
    }
    else
    {
        // Extrapolate

        if (extrap_method == Linear_Extrapolation)
        {
            return(linear_evaluate(xvalue));
        }    
    }
}














/**
 * \brief Initialize the necessary values for natural spline interpolation.
 * 
 * \details Compute the second derivatives of the interpolating 
 *          function at the tabulated points. By "natural" spline we mean
 *          that the second derivatives at the boundary values are assumed
 */

template <class Type>
void TabulatedFunction<Type>::init_spline()
{
    double p, qn, sig, un;

    // Set the lower boundary to be "natural", i.e. zero second derivative

    DDy[0] = temp[0] = 0.0;

    // Now the decomposition loop of the tridiagonal algorithm. DDy and u
    // are used for temporary storage of the decomposed factors.

    for (unsigned long i=1; i <= Nvalues-2; i++)
    { 
        sig = (x[i]-x[i-1]) / (x[i+1]-x[i-1]);
        p = sig * DDy[i-1] + 2.0;
        DDy[i] = (sig-1.0) / p;
        temp[i] =   (y[i+1]-y[i]) / (x[i+1]-x[i]) - (y[i]-y[i-1]) / (x[i]-x[i-1]);
        temp[i] = (6.0 * temp[i] / (x[i+1]-x[i-1]) - sig * temp[i-1]) / p;
    }

    // Also set the upper boundary to be the "natural" boundary

    qn = un = 0.0; 
    DDy[Nvalues-1]=(un - qn * temp[Nvalues-2]) / (qn * DDy[Nvalues-2] + 1.0);

    // Now the backsubstitution loop of the tridiagonal algorithm
    // Pitfall: don't declare k unsigned, otherwise it will never reach
    // the value -1 to stop the loop.

    for (long k = Nvalues-2; k >= 0; k--)
    { 
        DDy[k] = DDy[k] * DDy[k+1] + temp[k]; 
    }
}










/**
 * \details Computes the natural cubic spline interpolation, given
 *          the tabulated values of the function, and given the array
 *          of second derivatives computed by init_spline().
 *          "Natural" here means that the second derivatives of the function
 *          at the tabulated boundary values are assumed to be zero. 
 */

template <class Type>
double TabulatedFunction<Type>::spline_evaluate(const double xvalue)
{
    unsigned long klo, khi, k;
    double h,b,a,yvalue;

    // First find the right place in the table by means of bisection.

    klo = 0; khi = Nvalues-1;
    while (khi-klo > 1)
    { 
        k = (khi+klo) >> 1;
        if (x[k] > xvalue) 
        { 
            khi = k;
        }
        else 
        {
            klo = k;
        }
    }

    // Define a few handy abbriviations

    h = x[khi]-x[klo];
    a = (x[khi]-xvalue)/h;
    b = (xvalue-x[klo])/h;

    // Evaluate cubic spline

    yvalue = a*y[klo] + b*y[khi] + ((a*a*a-a)*DDy[klo]+(b*b*b-b)*DDy[khi]) * (h*h)/6.0;

    return(yvalue);
}













/**
 * \brief Perform linear inter- or extrapolation of a tabulated function
 */ 

template <class Type>
double TabulatedFunction<Type>::linear_evaluate(const double xvalue)
{
    long i,j;      // Indices of the two points, defining the linear relation
    long m;        // Index of middle point

    // Out of range, or on the upper border

    if (xvalue >= x[Nvalues-1])
    {
        i = Nvalues-2;
        j = Nvalues-1;
        return(y[i] + (y[j]-y[i])/(x[j]-x[i]) * (xvalue - x[i]));
    }

    // Out of range, or on the lower border

    if (xvalue <= x[0])
    {
        i = 0;
        j = 1;
        return(y[i] + (y[j]-y[i])/(x[j]-x[i]) * (xvalue - x[i]));
    }

    // In tabulated range. First locate neighbours by bisection.

    i = 0;                 // We already checked the lower and upper
    j = Nvalues-1;         // borders.

    while (j - i > 1)
    {
        m = (j + i) >> 1;          // middle point

        if (xvalue >= x[m])
        {
            i = m;
        }
        else
        {
            j = m;
        }
    }

  return(y[i] + (y[j]-y[i])/(x[j]-x[i]) * (xvalue - x[i]));

}













/**
 * \brief: Integrate the tabulated function. 
 *
 * \details The integrand is computed via inter/extrapolation.
 *
 * \param lower: Lower integration boundary
 * \param upper: Upper integration boundary
 *
 * \return  the integral from lower to upper.
 */

template <class Type>
double TabulatedFunction<Type>::integrate(double lower, double upper)
{
    const int JMAX = 20;
    double x,tnm,sum,del,olds;
    double s = 0.0;
    int it,j,k;

    // Init old s to a number that is unlikely to the average of the function
    // at its endpoints

    olds = -1.0e30;

    // Iterate with finer trapezium coverages until converged

    for (j = 1; j <= JMAX; j++)
    {
        if (j == 1)
        {
            s = 0.5 * (upper - lower);
        }
        else
        {
            for (it = 1, k = 1; k < j-1; k++) it <<= 1;
            tnm = it;
            del = (upper - lower)/tnm;
            x = lower + 0.5 * del;
            for (sum = 0.0, k=1; k <= it; k++, x+= del)
            {
                sum += x;
            }
            s = 0.5 * (s + (upper-lower)*sum / tnm);
        }

        // Avoid spurious early convergence, that's why the j > 5

        if (j > 5)
        {
            if (fabs(s-olds) < accuracy * fabs(olds) || (s == 0.0 && olds == 0.0))
            {
                return s;
            }
        }
      
        olds = s;
    }

    // If we get here, the convergence failed.

    Log.error("TabulatedFunction::Integrate()): no convergence.");
    exit(1);

}













/**
 * \brief  Set the fractional accuracy which determines when a convergence
 *         procedure (e.g. to integrate) has to stop.
 */

template <class Type>
void TabulatedFunction<Type>::setAccuracy(const double accur)
{
    if (accur > 0.0)
    {
        accuracy = accur;
    }
    else
    {
        Log.error("TabulatedFunction::setAccuracy(): accuracy <= 0");
        exit(1);
    }
}




#endif
