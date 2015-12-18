
#include "timeticker.h"



TimeTicker::TimeTicker()
{
    // Put the tick interval at the largest possible values.
    // This is a good value for classes inheriting from TimeTicker, but that 
    // do not care about a tick interval.

    tickInterval = numeric_limits<double>::max();
}






TimeTicker::~TimeTicker()
{
    
}





double TimeTicker::getTickInterval()
{
    return tickInterval;
}