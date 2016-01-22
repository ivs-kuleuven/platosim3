#ifndef STARCATALOG_H
#define STARCATALOG_H

#include <string>
#include <vector>
#include "logger.h"
#include "configurationparameters.h"

using namespace std;



class StarCatalog
{
    public:

        StarCatalog(ConfigurationParameters configurationParameters);
        ~StarCatalog();

        void getStarsWithinRadiusFrom(StarCatalog &starCatalog);
        void computeSkyBackground(double alpha, double delta);

        virtual void configureWithFile(string fileName);

    protected:


    private:

        string starCatalogFileName;
        string outputFileName;;                      // HDF5 file, including full path

        long Nstars;
        vector<long> starID;
        vector<double> rightAscension;
        vector<double> declination;
        vector<double> Vmag; 

};



#endif