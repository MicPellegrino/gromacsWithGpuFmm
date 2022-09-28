#ifndef _BK_FMM_local2local_derivative_hpp_
#define _BK_FMM_local2local_derivative_hpp_

#include "cuda_keywords.hpp"

namespace gmx_gpu_fmm{

template <typename CoefficientMatrix>
CUDA
void DxL2L_reference(
        const CoefficientMatrix & mu,
        CoefficientMatrix & dmu,
        const ssize_t p)
{
    typedef typename CoefficientMatrix::complex_type complex_type;
    typedef typename complex_type::value_type scalar_type;

    for (ssize_t l = 0; l < p; ++l)
    {
        for (ssize_t m = 0; m <= l; ++m)
        {
            dmu(l, m) = (mu.get_vectorized(l + 1, m + 1) - mu.get_vectorized(l + 1, m - 1)) * scalar_type(-0.5);
            //dmu(l, m) = (mu.get(l + 1, m + 1) - mu.get(l + 1, m - 1)) * scalar_type(-0.5);
        }

    }
    for (ssize_t m = 0; m <= p; ++m)
        dmu(p, m) = complex_type(0.);
}

template <typename CoefficientMatrix>
CUDA
void DyL2L_reference(
        const CoefficientMatrix & mu,
        CoefficientMatrix & dmu,
        const ssize_t p)
{
    typedef typename CoefficientMatrix::complex_type complex_type;
    typedef typename complex_type::value_type scalar_type;

    for (ssize_t l = 0; l < p; ++l) {
        for (ssize_t m = 0; m <= l; ++m) {
            complex_type tmp = (mu.get(l + 1, m + 1) + mu.get(l + 1, m - 1)) * scalar_type(0.5);
            dmu(l, m) = complex_type(-tmp.imag(), tmp.real());
        }
    }
    for (ssize_t m = 0; m <= p; ++m)
        dmu(p, m) = complex_type(0.);
}

template <typename CoefficientMatrix>
CUDA
void DzL2L_reference(
        const CoefficientMatrix & mu,
        CoefficientMatrix & dmu,
        const ssize_t p)
{
    typedef typename CoefficientMatrix::complex_type complex_type;

    for (ssize_t l = 0; l < p; ++l) {
        for (ssize_t m = 0; m <= l; ++m)
            dmu(l, m) = -mu.get(l + 1, m);
    }
    for (ssize_t m = 0; m <= p; ++m)
        dmu(p, m) = complex_type(0.);
}

}//namespace end

#endif
