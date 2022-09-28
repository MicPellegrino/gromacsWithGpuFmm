#include "timer.hpp"

namespace gmx_gpu_fmm{

Timer_::Timer_(){}

inline
double Timer_::timestamp()
{
    struct timeval tp;
    gettimeofday(&tp, nullptr);
    return double(tp.tv_sec) + tp.tv_usec / 1000000.;
}

void Timer_::start()
{
    start_time = timestamp();
}

void Timer_::stop()
{
    stop_time = timestamp();
}

double Timer_::get_time()
{
    return  stop_time - start_time;
}

}//namespace end
