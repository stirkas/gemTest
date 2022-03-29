#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>

#include "timer.hpp"

using namespace std;

void TimerClass::savetimes(string filename)
{
    ofstream timefile;
    timefile.open(filename);
    for(int i = 0; i < times.size(); ++i)
    {
        timefile << times[i].count() << endl;
    }
    timefile.close();
}


void TimerClass::startimer()
{
    t1 = Clock::now();
}

extern "C" void TimerClass::stoptimer()
{
    t2 = Clock::now();
    times.push_back(t2-t1);
}

