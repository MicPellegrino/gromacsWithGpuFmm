#include "managed.hpp"
#include "architecture.hpp"

namespace gmx_gpu_fmm{

template <typename arch>
void* Managed<arch>::operator new (size_t len)
{
    void *ptr;
    arch::custom_alloc(ptr,len);
    //printf("+++++++++++++++++++++++\n");
    return ptr;
}

template <typename arch>
void Managed<arch>::operator delete(void *ptr){
    arch::custom_free(ptr);
    //printf("---------------------\n");
}

template <typename arch>
void* Managed<arch>::operator new[] (size_t len)
{
    void *ptr;
    arch::custom_alloc(ptr,len);
    //printf("+++++++++++++++++++++++[]\n");
    return ptr;
}

template <typename arch>
void Managed<arch>::operator delete[](void *ptr){
    arch::custom_free(ptr);
    //printf("---------------------[]\n");
}

template class Managed<Device<double> >;
template class Managed<Device<float> >;
template class Managed<Host<double> >;
template class Managed<Host<float> >;



}//namespace end
