
#ifndef TABULATEDFUNCTION_H
#define TABULATEDFUNCTION_H

#include <vector>

using namespace std;


class TabulatedFunction
{
    public:

        TabulatedFunction(double *xValues, double *yValues, int Nvalues);
        ~TabulatedFunction();

        double integral(double lowerLimit, double upperLimit);
        double operator()(double xValue);
        vector<double> operator()(vector<double> &xValues);

        void setInterpolationMethod();
        void setIntegrationMethod();


    protected:

    private:

        vector<double> xvalues;
        vector<double> yvalues;
};



#endif
