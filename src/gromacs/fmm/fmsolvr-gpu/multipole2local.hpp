#ifndef _BK_FMM_multipole2local_hpp_
#define _BK_FMM_multipole2local_hpp_

#include "parity_sign.hpp"
#include "particle2local.hpp"
#include "bit_manipulator.hpp"
#include "architecture.hpp"

namespace gmx_gpu_fmm{

template <typename CoefficientMatrix, typename CoefficientMatrixDouble = CoefficientMatrix>
class M2L_Operator : public CoefficientMatrix
{
public:

    typedef typename CoefficientMatrix::architecture arch;
    typedef Bitit<32,5,arch> Bitset;
    Bitset* bitset;

    template <typename Real3>
    M2L_Operator( 
        const Real3 & xyz,
        const size_t p) : CoefficientMatrix(2 * p), cp(p)
    {
        tmp = new CoefficientMatrixDouble(2 * p);
        init(xyz,p);
    }

    template <typename Real3>
    M2L_Operator(
        const Real3 & expansion_point_source,
        const Real3 & expansion_point_target,
        const size_t p) : cp(p)

    {
        tmp = new CoefficientMatrixDouble(2 * p);
        init(expansion_point_source - expansion_point_target, p);
    }

    ~M2L_Operator()
    {
        //printf("pointer    --- >  %p\n", tmp);
        delete tmp;
    }

    void init_bitset(Bitset* bitset_mem)
    {
        bitset = bitset_mem;
    }

    void set_bitset()
    {
        /*
        for(ssize_t l = 0; l <= 2*cp; ++l)
        {
            for(ssize_t m = 0; m <= l; ++m)
            {
                printf("(%e %e) ", this->get(l,m).real(), this->get(l,m).imag());
            }
            printf("\n");
        }
        */
        size_t index_global = 0;
        for(ssize_t l = 0; l <= cp; ++l)
        {
            for(ssize_t m = 0; m <= l; ++m)
            {
                size_t index_local = 0;
                //printf("l %lu m %lu\n",l,m);
                //bitset[index_global].dump();
                for(ssize_t j = 0; j <= cp; ++j)
                {
                    for(ssize_t k = -j; k <= j; ++k)
                    {
                        //printf("(%e %e) ", this->get(l+j,m+k).real(), this->get(l+j,m+k).imag());

                        if(signbit(this->get(l+j,m+k).real()))
                        {
                            bitset[index_global].set(index_local);
                        }
                        ++index_local;

                        if(signbit(this->get(l+j,m+k).imag()))
                        {
                            bitset[index_global].set(index_local);
                        }
                        ++index_local;
                    }
                    //printf("\n");
                }
                //bitset[index_global].dump();
                //printf("\n");
                ++index_global;
            }
        }
    }

private:
    CoefficientMatrixDouble* tmp;

    template <typename Real3>
    void init(
        const Real3 & xyz,
        const size_t p)
    {
        typedef typename Real3::value_type Real;

        P2L(make_xyzq(xyz, double(1.0)), *tmp, 2 * p);

        this->recast(*tmp);
        this->populate_lower();
    }
    const size_t cp;
};


template <typename VCoefficientMatrix>
class M2L_Operators : public VCoefficientMatrix
{
public:

    M2L_Operators(){}

    M2L_Operators( const size_t p,
                   const size_t n) : VCoefficientMatrix(2*p, n)
    {}
};


// problematic to get the shifting vector xyz in the right orientation
template <typename Real3, typename CoefficientMatrix>
void M2L_reference_DONOTUSE(
        const CoefficientMatrix & omega,
        const Real3 & xyz,
        CoefficientMatrix & mu,
        const ssize_t p)
{
    M2L_Operator<CoefficientMatrix> B(xyz, p);
    M2L_reference(omega, B, mu, p);
}


// ensure the correct orientation of the shifting vector
template <typename Real3, typename CoefficientMatrix>
void M2L_reference(
        const CoefficientMatrix & omega_source,
        const Real3 & expansion_point_source,
        CoefficientMatrix & mu_target,
        const Real3 & expansion_point_target,
        const ssize_t p)
{
    M2L_Operator<CoefficientMatrix> B(expansion_point_source, expansion_point_target, p);
    M2L_reference(omega_source, B, mu_target, p);
}


// p^4, k=-j..j, with custom B operator
// classical M2L operation on actual *charges*, does include the (-1)^j factor
template <typename CoefficientMatrix, typename OperatorMatrix>
void M2L_reference(
        const CoefficientMatrix & omega,
        const OperatorMatrix & B,
        CoefficientMatrix & mu,
        const ssize_t p)
{
    typedef typename CoefficientMatrix::complex_type complex_type;

    for (ssize_t l = 0; l <= p; ++l)
    {
        for (ssize_t m = 0; m <= l; ++m)
        {
            complex_type mu_l_m(0.);
            for (ssize_t j = 0; j <= p; ++j)
            {
                complex_type mu_l_m_j(0.);
                for (ssize_t k = -j; k <= j; ++k)
                {
                    mu_l_m_j += B.get(j + l, k + m) * omega.get(j, k);
                }
                mu_l_m += toggle_sign_if_odd(j, mu_l_m_j);
            }

            mu(l, m) += toggle_sign_if_odd(0, mu_l_m);
        }
    }
}


// p^4, k=-j..j, with custom B operator
// raw == for operator construction (chargeless), does *not* include (-1)^{j,l} factor
template <typename CoefficientMatrix, typename OperatorMatrix>
void rawM2L_reference(
        const CoefficientMatrix & omega,
        const OperatorMatrix & B,
        CoefficientMatrix & mu,
        const ssize_t p)
{
    typedef typename CoefficientMatrix::complex_type complex_type;
    const ssize_t p_omega = std::min<ssize_t>(p, omega.p());
    const ssize_t p_B = B.p();
    const ssize_t p_mu = std::min<ssize_t>(p, mu.p());

    for (ssize_t l = 0; l <= p_mu; ++l) {
        for (ssize_t m = 0; m <= l; ++m) {
            complex_type mu_l_m(0.);
            const ssize_t j_max = std::min(p_omega, p_B - l);
            for (ssize_t j = 0; j <= j_max; ++j) {
                for (ssize_t k = -j; k <= j; ++k)
                    mu_l_m += B.get(j + l, k + m) * omega.get(j, k);
            }
            mu(l, m) += mu_l_m;
        }
    }
}


// classical M2L operation on actual *charges*, does include the (-1)^l factor
template <typename CoefficientMatrix>
void
__attribute__ ((noinline))
__attribute ((flatten))
    single_M2L(
        const CoefficientMatrix & omega,
        const CoefficientMatrix & B,
        CoefficientMatrix & mu,
        const ssize_t p)
{
    M2L_reference(omega, B, mu, p);
}


// p^4, k=-j..j, with custom B operator
// raw == for operator construction (chargeless), does *not* include (-1)^l factor
template <typename CoefficientMatrix>
void
__attribute__ ((noinline))
__attribute ((flatten))
    single_rawM2L(
        const CoefficientMatrix & omega,
        const CoefficientMatrix & B,
        CoefficientMatrix & mu,
        const ssize_t p)
{
    rawM2L_reference(omega, B, mu, p);
}

}//namespace end


#endif
