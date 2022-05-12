#pragma once

#include "timer.hpp"

#include <iostream>

using namespace std;

extern "C"
{
    void initptr(TimerClass *ptrtimer);


    void savetimesf(TimerClass *ptrtimer);


    void starttimerf(TimerClass *ptrtimer);


    void stoptimerf(TimerClass *ptrtimer);

    void testcase_();
 
}

