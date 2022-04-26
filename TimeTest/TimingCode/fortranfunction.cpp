#include "fortranfunction.h"

#include <iostream>

using namespace std;


void initptr(TimerClass *ptrtimer)
{
    if(NULL != ptrtimer)
    {
        TimerClass timer;
        ptrtimer = &timer;
    }

    
}
void savetimesf(TimerClass *ptrtimer)
{
    if(NULL != ptrtimer)
    {
        ptrtimer->savetimes("fortrantimes.txt");
    }
    
}
void starttimerf(TimerClass *ptrtimer)
{
    if(NULL != ptrtimer)
    {
        ptrtimer->startimer();
    }
    
}
void stoptimerf(TimerClass *ptrtimer)
{
    if(NULL!= ptrtimer)
    {
       ptrtimer->stoptimer(); 
    }
    
}