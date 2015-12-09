#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <string>
#include "detector.h"
#include "camera.h"
#include "telescope.h"
#include "platform.h"
#include "sky.h"
#include "configurationparameters.h"


using namespace std;


class Controller
{
    public:

        Controller(string inputFileName);
        ~Controller();
        virtual void run(double startingTime = 0.0);

    protected:

        
    private:

        double currentTime;

        Detector  *detector;
        Camera    *camera;
        Telescope *telescope;
        Platform  *platform;
        Sky       *sky;

};



#endif
