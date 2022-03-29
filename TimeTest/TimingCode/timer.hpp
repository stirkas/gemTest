#include <chrono>
#include <string>
#include <vector>

class TimerClass
{
    typedef std::chrono::high_resolution_clock Clock;

    public:
        TimerClass(){}
        ~TimerClass(){}
        void savetimes(std::string filename);
        void startimer();
        void stoptimer();
    private:
        std::vector<Clock::duration> times;
        Clock::time_point t1;
        Clock::time_point t2; 
};