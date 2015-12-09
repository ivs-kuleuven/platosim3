#ifndef JITTER_H
#define JITTER_H

#include <string>
#include <vector>



using namespace std;



class JitterGenerator
{
    public:

        JitterGenerator();
        ~JitterGenerator();

        void getNextYawPitchRoll(double &yaw, double &pitch, double &roll, double timeInterval);

    protected:


    private:

 };



#endif