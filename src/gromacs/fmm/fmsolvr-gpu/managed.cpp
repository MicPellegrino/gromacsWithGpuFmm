#include "managed.hpp"
#include "architecture.hpp"

namespace gmx_gpu_fmm{

template <typename arch>
void* Managed::operator new (size_t len)
{
    void *ptr;
    arch::custom_alloc(ptr,len);
    return ptr;
}

template <typename arch>
void Managed::operator delete(void *ptr){
    arch::custom_free(ptr);
}


template <typename arch>
void* Managed::operator new[] (size_t len)
{
    void *ptr;
    arch::custom_alloc(ptr,len);
    return ptr;
}

template <typename arch>
void Managed::operator delete[](void *ptr){
    arch::custom_free(ptr);
}


}//namespace end
