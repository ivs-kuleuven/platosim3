
#ifndef TIMETICKER_H
#define TIMETICKER_H

#include <limits>

using namespace std;



class TimeTicker
{
    
    public:

        TimeTicker();
        ~TimeTicker();

        double getTickInterval();

    protected:

        double tickInterval;       // [seconds]

    private:

    
};

#endif
