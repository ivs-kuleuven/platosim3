#ifndef ASYMMETRICALPOINTSPREADFUNCTION_H
#define ASYMMETRICALPOINTSPREADFUNCTION_H

#include <string>
#include <vector>
#include <algorithm>

// #include "armadillo"

#include "ConfigurationParameters.h"
#include "PointSpreadFunction.h"

using namespace std;





class AsymmetricalPointSpreadFunction : public PointSpreadFunction
{
    public:

        AsymmetricalPointSpreadFunction(ConfigurationParameters &configParam, HDF5File &hdf5File);

        void select(double xFP, double yFP);

        virtual void rotate(double angle);

    protected:

        virtual void configure(ConfigurationParameters &configParam);
};

#endif
