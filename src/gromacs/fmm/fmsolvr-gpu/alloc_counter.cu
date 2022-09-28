#include "alloc_counter.hpp"

namespace gmx_gpu_fmm{

int alloc_counter::allocs = 0;

void alloc_counter::up()
{
    alloc_counter::allocs++;
}

void alloc_counter::down()
{
    alloc_counter::allocs--;
}

int alloc_counter::get()
{
    return alloc_counter::allocs;
}

}//namespace end

