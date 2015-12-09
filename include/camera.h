#ifndef CAMERA_H
#define CAMERA_H

#include <string>
#include "telescope.h"
#include "subfield.h"



using namespace std;



class Camera
{
    public:

        Camera(ConfigurationParameters configurationParameters);
        ~Camera();

        void initPsf(SubField subField);
        void exposeSubField(SubField subfield);

    protected:


    private:

        Telescope telescope;

        double plateScale;             // [arcsec/mm]
        double internalTime;           // [s]

};



#endif
