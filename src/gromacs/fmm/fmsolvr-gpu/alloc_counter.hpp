#ifndef ALLOC_COUNTER_HPP
#define ALLOC_COUNTER_HPP

namespace gmx_gpu_fmm{

/*! \brief Used for counting the number of allocations
 *
 */
class alloc_counter
{
public:
    static int allocs;

    static void up();

    static void down();

    static int get();
};

}//namespace end
#endif // ALLOC_COUNTER_HPP
