
#include "Heartbeat.h"



Heartbeat::Heartbeat()
{
    // Put the tick interval at the largest possible values.
    // This is a good value for classes inheriting from Heartbeat, but that 
    // do not care about a tick interval.

    heartbeatInterval = numeric_limits<double>::max();
}






Heartbeat::~Heartbeat()
{
    
}





double Heartbeat::getHeartbeatInterval()
{
    return heartbeatInterval;
}
