#ifndef _BK_FMM_multipole2multipole_derivative_hpp_
#define _BK_FMM_multipole2multipole_derivative_hpp_

namespace gmx_gpu_fmm{

template <typename CoefficientMatrix>
void DxM2M_reference(
        const CoefficientMatrix & omega,
        CoefficientMatrix & domega,
        const ssize_t p)
{
    typedef typename CoefficientMatrix::complex_type complex_type;
    typedef typename complex_type::value_type scalar_type;

    if (p >= 0)
        domega(0, 0) = complex_type(0.);
    if (p >= 1)
        domega(1, 0) = complex_type(0.);
    for (ssize_t l = 1; l <= p; ++l) {
        for (ssize_t m = 0; m <= l - 2; ++m)
            domega(l, m) = (omega.get(l - 1, m - 1) - omega.get(l - 1, m + 1)) * scalar_type(0.5);
        for (ssize_t m = std::max<ssize_t>(l - 1, 1); m <= l; ++m)
            domega(l, m) = omega.get(l - 1, m - 1) * scalar_type(0.5);
    }
}

template <typename CoefficientMatrix>
void DyM2M_reference(
        const CoefficientMatrix & omega,
        CoefficientMatrix & domega,
        const ssize_t p)
{
    typedef typename CoefficientMatrix::complex_type complex_type;
    typedef typename complex_type::value_type scalar_type;

    if (p >= 0)
        domega(0, 0) = complex_type(0.);
    if (p >= 1)
        domega(1, 0) = complex_type(0.);
    for (ssize_t l = 1; l <= p; ++l) {
        for (ssize_t m = 0; m <= l - 2; ++m) {
            complex_type tmp = (omega.get(l - 1, m - 1) + omega.get(l - 1, m + 1)) * scalar_type(0.5);
            domega(l, m) = complex_type(tmp.imag(), -tmp.real());
        }
        for (ssize_t m = std::max<ssize_t>(l - 1, 1); m <= l; ++m) {
            complex_type tmp = omega.get(l - 1, m - 1) * scalar_type(0.5);
            domega(l, m) = complex_type(tmp.imag(), -tmp.real());
        }
    }
}

template <typename CoefficientMatrix>
void DzM2M_reference(
        const CoefficientMatrix & omega,
        CoefficientMatrix & domega,
        const ssize_t p)
{
    typedef typename CoefficientMatrix::complex_type complex_type;

    if (p >= 0)
        domega(0, 0) = complex_type(0.);
    for (ssize_t l = 1; l <= p; ++l) {
        for (ssize_t m = 0; m < l; ++m)
            domega(l, m) = omega.get(l - 1, m);
        domega(l, l) = complex_type(0.);
    }
}

}//namespace end

#endif
