#ifndef _BK_FMM_local2local_hpp_
#define _BK_FMM_local2local_hpp_

#include <complex>
#include "particle2multipole.hpp"

namespace gmx_gpu_fmm{

template <typename CoefficientMatrix, typename CoefficientMatrixDouble = CoefficientMatrix>
class L2L_Operator : public CoefficientMatrix
{
public:
    template <typename Real3>
    L2L_Operator(
        const Real3 & xyz,
        const size_t p) : CoefficientMatrix(p)
    {
        init(xyz,p);
    }

    template <typename Real3>
    L2L_Operator(
        const Real3 & expansion_point_source,
        const Real3 & expansion_point_target,
        const size_t p): CoefficientMatrix(p)
    {
        init(expansion_point_target - expansion_point_source, p);
    }

    template <typename Real3>
    void init_new(
        const Real3 & xyz,
        const size_t p)
    {
        typedef typename Real3::value_type Real;
        (*static_cast<CoefficientMatrix *>(this)).zero();
        P2M_nodivtable(Real(1.0), xyz, *static_cast<CoefficientMatrix *>(this), p);
        this->populate_lower();
    }

    private:
    template <typename Real3>
    void init(
        const Real3 & xyz,
        const size_t p)
    {
        typedef typename Real3::value_type Real;
        P2M_nodivtable(Real(1.0), xyz, *static_cast<CoefficientMatrix *>(this), p);
        this->populate_lower();
    }


};


// ensure the correct orientation of the shifting vector
template <typename Real3, typename CoefficientMatrix>
void L2L_reference(
        const CoefficientMatrix & mu_source,
        const Real3 & expansion_point_source,
        CoefficientMatrix & mu_target,
        const Real3 & expansion_point_target,
        const ssize_t p,
        const ssize_t p_out)
{
    assert(p >= p_out);
    L2L_Operator<CoefficientMatrix> C(expansion_point_source, expansion_point_target, p);
    L2L_reference(mu_source, C, mu_target, p, p_out);
}


// p^4, k=-j..j
template <typename Real3, typename CoefficientMatrix>
void L2L_reference(
        const CoefficientMatrix & mu_in,
        const Real3 & xyz,
        CoefficientMatrix & mu_out,
        const ssize_t p)
{
    L2L_Operator<CoefficientMatrix> C(xyz, p);
    L2L_reference(mu_in, C, mu_out, p, p);
}


// p_in^2*p_out^2, k=-j..j
template <typename Real3, typename CoefficientMatrix>
void L2L_reference(
        const CoefficientMatrix & mu_in,
        const Real3 & xyz,
        CoefficientMatrix & mu_out,
        const ssize_t p,
        const ssize_t p_out)
{
    assert(p >= p_out);
    L2L_Operator<CoefficientMatrix> C(xyz, p);
    L2L_reference(mu_in, C, mu_out, p, p_out);
}


// p^4 (actually p^2 * p_out^2), k=-j..j
template <typename CoefficientMatrix>
void L2L_reference(
        const CoefficientMatrix & mu_in,
        const L2L_Operator<CoefficientMatrix> & C,
        CoefficientMatrix & mu_out,
        const ssize_t p,
        const ssize_t p_out_)
{
    typedef typename CoefficientMatrix::complex_type complex_type;
    // fail if mu_out(0,0) cannot be computed correctly due to insufficient input order
    //printf("before assert %lu   %lu  %lu \n",p, mu_in.p(), C.p());
    assert(p <= std::min<ssize_t>(mu_in.p(), C.p()));
    // limit output order to values that can be stored and computed correctly
    const ssize_t p_out = std::min<ssize_t>(std::min(p_out_, p), mu_out.p());

    for (ssize_t l = 0; l <= p_out; ++l) {
        for (ssize_t m = 0; m <= l; ++m) {
            complex_type mu_l_m(0.);
            for (ssize_t j = l; j <= p; ++j) {
                const ssize_t k_min = m - (j - l);
                const ssize_t k_max = m + (j - l);
                for (ssize_t k = k_min; k <= k_max; ++k)
                    mu_l_m += C.get(j - l, k - m) * mu_in.get(j, k);
            }
            mu_out(l, m) += mu_l_m;
        }
    }
}


// ensure the correct orientation of the shifting vector
template <typename Real3, typename CoefficientMatrix>
void L2L_reference(
        const CoefficientMatrix & mu_source,
        const Real3 & expansion_point_source,
        CoefficientMatrix & mu_target,
        const Real3 & expansion_point_target,
        const ssize_t p)
{
    L2L_Operator<CoefficientMatrix> C(expansion_point_source, expansion_point_target, p);
    L2L_reference(mu_source, C, mu_target, p, p);
}

}//namespace end

#endif


