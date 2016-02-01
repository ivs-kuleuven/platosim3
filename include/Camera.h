#ifndef CAMERA_H
#define CAMERA_H

#include <string>

#include "Logger.h"
#include "TimeTicker.h"
#include "HDF5Writer.h"
#include "Telescope.h"


using namespace std;



class Camera : public TimeTicker, HDF5Writer
{
    public:

        Camera(ConfigurationParameters configurationParameters);
        ~Camera();

        void initPsf(SubField subField);
        void exposeSubField(Dectector &detector);

    protected:


    private:

        Telescope telescope;

        double plateScale;             // [arcsec/mm]
        double internalTime;           // [s]

};



#endif
