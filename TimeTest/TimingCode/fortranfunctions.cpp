#include "timer.hpp"

#include <iostream>

using namespace std;

extern "C"
{
    void initptr()
    {
        TimerClass timer;

        TimerClass *ptrtimer;

        ptrtimer = &timer;
    }

    void savetimesf(TimerClass *ptrtimer)
    {
        ptrtimer->savetimes("fortrantimes.txt");
    }

    void starttimerf(TimerClass *ptrtimer)
    {
        ptrtimer->startimer();
    }

    void stoptimerf(TimerClass *ptrtimer)
    {
        ptrtimer->stoptimer();
    }
}

