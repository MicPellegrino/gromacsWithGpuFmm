#include "global_functions.hpp"

namespace gmx_gpu_fmm{

int get_env_int(int dflt, const char* vrnm)
{
    if (const char * val = getenv(vrnm))
        return atoi(val);
    else
        return dflt;
}

size_t boxes_above_depth(size_t d)
{
    // sum(0, d-1, 8^i)
    return ((size_t(1) << (3 * d)) - 1 ) / 7;
}

size_t boxes_on_depth(size_t d)
{
    // 8^d
    return size_t(1) << (3 * d);
}

size_t make_boxid(size_t x, size_t y, size_t z, unsigned depth)
{
    return (z << (depth * 2)) + (y << depth) + x;
}

}//namespace end
