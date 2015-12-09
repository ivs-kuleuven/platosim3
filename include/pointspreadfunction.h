#ifndef PSF_H
#define PSF_H

#include <string>
#include <vector>
#include "configurationparameters.h"

using namespace std;


class PointSpreadFunction
{
    public:

        PointSpreadFunction(ConfigurationParameters configurationParameters);
        ~PointSpreadFunction();


    protected:


    private:

        void select();
        void rotate();
        void rebin();

};



#endif