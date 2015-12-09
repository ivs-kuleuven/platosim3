#ifndef SUBFIELD_H
#define SUBFIELD_H

#include <string>
#include <vector>

#include "configurationparameters.h"
#include "starcatalog.h"


using namespace std;


class SubField
{
    public:

        SubField();
        ~SubField();

        double getDistanceFromOpticalAxisToFieldCenter();
        void addFlux(double xCoord, double yCoord, double flux);
        void convolveWithPsf(double **psf);
        void multiply(double **array);
        void rebin();
        void reset();


    protected:


    private:


};



#endif