#include "timer.hpp"

#include <iostream>

using namespace std;

int main()
{
    int counter = 0;

    TimerClass timer;
    for(int j = 0; j < 101; ++j)
    {
        timer.startimer();
        for(int i = 0; i <j; ++i)
        {
            counter ++;
        }

        timer.stoptimer();
        counter = 0;
    }

    timer.savetimes("testfile.txt");

    return 0;
}