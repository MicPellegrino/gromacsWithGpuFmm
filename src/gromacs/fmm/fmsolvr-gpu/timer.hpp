#ifndef _BK_TIMER_HPP_
#define _BK_TIMER_HPP_

#include <ctime>
#include <sys/time.h>
#include <map>
#include <vector>
#include <algorithm>


namespace gmx_gpu_fmm{

struct Timer_
{
    double start_time;
    double stop_time;

    Timer_();

    double timestamp();

    void start();

    void stop();

    double get_time();
};

class Exe_time{

public:
    double _walltime;
    double steps, actual_step;
    double offset;

    Exe_time():_walltime(0.0), offset(-2.0), actual_step(-1.0)
    {

    }

    void init(int nsteps)
    {
        actual_step = 0.0;
        _walltime = 0.0;
        steps = nsteps;
        for(int i = 0; i < std::min((int)steps, 22); i++)
        {
            offset = (double)i - 1.0;
        }
        //printf("offset %f\n", offset);

    }

    void init(std::string n)
    {
        _name[n] = n;
        _time[n] = 0.0;
    }

    void set(std::string n, double t)
    {
        _name[n] = n;
        if(actual_step > offset || actual_step == -1.0)
        {
            _time[n] += t;
            _walltime += t;
        }

    }

    std::string name(std::string n){
        return _name[n];
    }

    double time(std::string n){

        return _time[n]/(actual_step - offset);
    }

    double walltime(){
        return _walltime/(actual_step - offset);
    }

    void add_step()
    {
        actual_step++;
    }

    void start(std::string n)
    {
#ifdef CUDADEBUG
        timers[n].start();
#endif
    }
    void stop(std::string n)
    {
#ifdef CUDADEBUG
        timers[n].stop();
        set(n,timers[n].get_time());
#endif
    }

private:
    std::map<std::string, std::string> _name;
    std::map<std::string, double> _time;
    //std::map<std::string, std::vector<double>> _times;
    std::map<std::string, Timer_> timers;
};

}//namespace end

#endif
