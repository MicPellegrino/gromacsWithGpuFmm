#ifndef _BK_FMM_fmm_complex_hpp_
#define _BK_FMM_fmm_complex_hpp_

#include "cuda_keywords.hpp"
#include "cuda_complex.hpp"

namespace gmx_gpu_fmm{

template <typename T>
struct COMPLEX_GENERATOR
{
    typedef complex<T> type;
    typedef type* pointer;
};

}//namespace end

#endif

/*
namespace fmm
{

    template <typename Tp>

#ifdef __CUDACC__
    class complex : public cuda::complex<Tp>
#else
    class complex : public std::complex<Tp>
#endif

    {
        
#ifdef __CUDACC__
    typedef cuda::complex<Tp> Base;
#else
    typedef std::complex<Tp> Base;
#endif
        
        
    public:

        CUDA
        complex() : Base() { }
         
        CUDA
        complex(const Tp & r) : Base(r) {}
        
        CUDA
        complex(const Tp & r, const Tp & i) : Base(r, i) { }
        
        CUDA
        complex(const Base & c) : Base(c) { }
        
        CUDA
        complex<Tp> & operator *= (const complex<Tp>& z)
        {
            
            const Tp re = this->real() * z.real() - this->imag() * z.imag();
            const Tp im = this->real() * z.imag() + this->imag() * z.real();
            this->real(re);
            this->imag(im);
            return *this;
        }


    };
}
*/
