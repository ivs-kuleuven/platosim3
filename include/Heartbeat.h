
#ifndef HEARTBEAT_H
#define HEARTBEAT_H

#include <limits>

using namespace std;



class Heartbeat
{
    
    public:

        Heartbeat();
        ~Heartbeat();

        virtual double getHeartbeatInterval();

    protected:

        double heartbeatInterval;       // [seconds]

    private:

    
};

#endif
