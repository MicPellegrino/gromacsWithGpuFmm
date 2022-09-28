#ifndef _BK_FMM_multipole2multipole_hpp_
#define _BK_FMM_multipole2multipole_hpp_

#ifdef __CUDACC__
#include "cuda_complex.hpp"
#else
#include <complex>
#endif

#include "particle2multipole.hpp"
#include "xyzq.hpp"

namespace gmx_gpu_fmm{

template <typename CoefficientMatrix,  typename CoefficientMatrixDouble = CoefficientMatrix>
class M2M_Operator : public CoefficientMatrix
{
public:
    template <typename Real3>
    M2M_Operator(
        const Real3 & xyz,
        const size_t p) : CoefficientMatrix(p)
    {
        //printf("M2M\n");
        init(xyz,p);
    }

    template <typename Real3>
    M2M_Operator(
        const Real3 & expansion_point_source,
        const Real3 & expansion_point_target,
        const size_t p)
    : CoefficientMatrix(p)
    {
        tmp = new CoefficientMatrixDouble(p);
        init(expansion_point_source - expansion_point_target, p);
    }

    ~M2M_Operator()
    {
        delete tmp;
    }

private:
    CoefficientMatrixDouble* tmp;

    template<typename Real3>
    void init (const Real3 & xyz,
               const size_t p)
    {
        typedef typename Real3::value_type Real;
        P2M(make_xyzq(xyz, double(1.0)), *tmp, p);

        this->recast(*tmp);
        //dump(*temp,p);
        //this->populate_lower();
    }
};


// p^4
template <typename Real3, typename CoefficientMatrix>
void M2M_reference(
        const CoefficientMatrix & omega_source,
        const Real3 & expansion_point_source,
        CoefficientMatrix & omega_target,
        const Real3 & expansion_point_target,
        const ssize_t p)
{
    M2M_Operator<CoefficientMatrix> A(expansion_point_source, expansion_point_target, p);
    M2M_reference(omega_source, A, omega_target, p);
}


// p^4, deprecated interface (xyz may easily be misoriented)
template <typename Real3, typename CoefficientMatrix>
__attribute__ ((deprecated))
void M2M_reference(
        const CoefficientMatrix & omega_in,
        const Real3 & xyz,
        CoefficientMatrix & omega_out,
        const ssize_t p)
{
    M2M_Operator<CoefficientMatrix> A(xyz, p);
    M2M_reference(omega_in, A, omega_out, p);
}


// p^4, k=-j..j
// this function may also be used to perform P2M of (mono-/)di-/quadrupoles ...
// by feeding them in via omega_in with omega_in.p() \in {1,2,3,...}
template <typename CoefficientMatrix>
void M2M_reference(
        const CoefficientMatrix & omega_in,
        const M2M_Operator<CoefficientMatrix> & A,
        CoefficientMatrix & omega_out,
        const ssize_t p)
{
    typedef typename CoefficientMatrix::complex_type complex_type;
    const ssize_t p_in = omega_in.p();
    const ssize_t p_out = std::min<ssize_t>(p, std::min(omega_out.p(), A.p()));

    for (ssize_t l = 0; l <= p_out; ++l) {
        const ssize_t j_max = std::min(l, p_in);
        for (ssize_t m = 0; m <= l; ++m) {
            complex_type omega_l_m(0.);
            for (ssize_t j = 0; j <= j_max; ++j) {
                const ssize_t k_min = std::max(-j, m - (l - j));
                const ssize_t k_max = std::min(+j, m + (l - j));
                //printf("l %lu m %lu j %lu !!!!! j_max %lu+++++++++ k_min %lu k_max %lu \n",(int)l,(int)m,(int)j,(int)j_max,  (int)k_min,(int)k_max);

                //printf("%lu%lu%lu %lu%lu%lu\n",(int)l,(int)m,(int)j,(int)j_max,  (int)k_min,(int)k_max);
                for (ssize_t k = k_min; k <= k_max; ++k)
                    omega_l_m += A.get(l - j, m - k) * omega_in.get(j, k);
            }
            omega_out(l, m) += omega_l_m;
        }
    }
}

}//namespace end


#endif
