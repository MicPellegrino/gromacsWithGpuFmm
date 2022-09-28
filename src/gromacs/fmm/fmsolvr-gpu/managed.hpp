#ifndef MANAGED_HPP
#define MANAGED_HPP

namespace gmx_gpu_fmm{

template<typename arch>
class Managed{

public:

    void *operator new (size_t len);
    void operator delete(void *ptr);

    void *operator new[] (size_t len);
    void operator delete[](void *ptr);
};

}//namespace end

#endif // MANAGED_HPP
