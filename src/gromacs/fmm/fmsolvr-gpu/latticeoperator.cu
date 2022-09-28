#ifndef _BK_FMM_latticeoperator_hpp_
#define _BK_FMM_latticeoperator_hpp_

#include <limits>
#include "latticeoperator.hpp"
#include "lstar_operator.hpp"
#include "ostar_operator.hpp"

namespace gmx_gpu_fmm{

// S_L operator *** OVERWRITES OUTPUT ***
template <typename CoefficientMatrix>
void ScaleL_operator(
        const CoefficientMatrix & mu_in,
        const size_t ws,
        CoefficientMatrix & mu_out,
        const size_t p)
{
    typedef typename CoefficientMatrix::complex_type complex_type;
    typedef typename complex_type::value_type Real;
    Real third = reciprocal<Real>(2 * ws + 1);
    Real scale = third;
    for (ssize_t l = 0; l <= ssize_t(p); ++l) {
        for (ssize_t m = 0; m <= l; ++m) {
            mu_out(l, m) = scale * mu_in(l, m);
        }
        scale *= third;
    }
}

template <typename CoefficientMatrix>
void CopyTriangle(
        const CoefficientMatrix & in,
        CoefficientMatrix & out,
        const size_t p)
{
    for (ssize_t l = 0; l <= ssize_t(p); ++l) {
        for (ssize_t m = 0; m <= l; ++m) {
            out(l, m) = in(l, m);
        }
    }
}

// lattice contribution from the box centered at (x,y,z)
template <typename Real, typename CoefficientMatrix>
inline void Lattice(
        const Real x,
        const Real y,
        const Real z,
        CoefficientMatrix & lattice,
        const size_t p)
{
    typedef Real real_type;
    typedef typename CoefficientMatrix::complex_type complex_type;

    real_type dist_squared = x * x + y * y + z * z;
    real_type recip_dist = rsqrt(dist_squared);
    real_type recip_dist_squared = reciprocal(dist_squared);
    real_type recip_dist_to_the_4 = recip_dist_squared * recip_dist_squared;

    // lattice_0_0  (for the first iteration)
    complex_type lattice_mplus2_mplus2(recip_dist);

    // lattice_0_0 upto lattice_p_p-1
    for (size_t m = 0; m < p; m += 2) {
        // lattice_m_m  (from previous iteration)
        complex_type lattice_m_m = lattice_mplus2_mplus2;
        lattice(m, m) += lattice_m_m;

        // lattice_m+2_m+2  (for the next iteration)
        lattice_mplus2_mplus2 = real_type((2 * m - 1) * (2 * m - 3)) * recip_dist_to_the_4 * complex_type(x, y) * lattice_m_m;

        // lattice_m+2_m upto lattice_p_m
        complex_type lattice_lminus2_m = lattice_m_m;
        for (size_t l = m + 2; l <= p; l += 2) {
            // lattice_l_m
            complex_type lattice_l_m = -real_type((l - m - 1) * (l + m - 1)) * recip_dist_squared * lattice_lminus2_m;
            lattice(l, m) += lattice_l_m;
            lattice_lminus2_m = lattice_l_m;
            if (l == 2 && m == 0)
                std::cout << x << " " << y << " " << z << " " << lattice_l_m << std::endl;
        }
    }
    // lattice_p_p
    lattice(p, p) += lattice_mplus2_mplus2;
}

// compute lattice operator
// (input lattice must be either zeroed or contain proper coefficients to be added to)
template <typename Real3, typename CoefficientMatrix>
void brokenLatticeOperator(
        Real3 & a,
        Real3 & b,
        Real3 & c,
        CoefficientMatrix & lattice,
        const size_t p)
{
    // first shell
    // ws == 1: 27 * 26 contributers
    const int outer = 4;
    const int inner = 1;
    int count = 0;
    for (int i = -outer; i <= +outer; ++i)
        for (int j = -outer; j <= +outer; ++j)
            for (int k = -outer; k <= +outer; ++k)
                if (i < -inner || i > inner || j < -inner || j > inner || k < -inner || k > inner) {
                    ++count;
                    Real3 m = a * i + b * j + c * k;
//                    std::cout << m.x << " " << m.y << " " << m.z << std::endl;
                    Lattice(m.x, m.y, m.z, lattice, p);
                }
//    std::cout << "Lattice eval: " << count << std::endl;
}

template <typename CoefficientMatrix>
size_t __attribute__((noinline)) LatticeOperatorIterate(
        const CoefficientMatrix & O_star,
        const CoefficientMatrix & L_star,
        const ssize_t ws,
        CoefficientMatrix & L_i,
        const size_t p,
        const size_t maxlayers,
        ssize_t check_convergence_at,
        const double maxerror = std::numeric_limits<double>::min())
{
    typedef typename CoefficientMatrix::value_type::value_type Real;
    Real oldval = L_star(check_convergence_at, 0).real();
    CoefficientMatrix L_pred_scaled(p);
    for (size_t i = 0; i < maxlayers; ++i) {
        if (i == 0)
            ScaleL_operator(L_star, ws, L_pred_scaled, p);
        else
            ScaleL_operator(L_i, ws, L_pred_scaled, p);
        L_pred_scaled.populate_lower();
        // compute L_{i+1}
        CopyTriangle(L_star, L_i, p);
        single_rawM2L(O_star, L_pred_scaled, L_i, p);
        // TODO: M2L operator with customizable loop stepping
        //       - no need to compute elements that are provable 0.0
        Real newval = L_i(check_convergence_at, 0).real();
        if (std::abs(oldval - newval) < maxerror)
            return i + 1;
        oldval = newval;
    }
    return maxlayers + 1;
}

// compute lattice operator
// 'lattice' will be overwritten, not accumulated into
template <bool enable_x, bool enable_y, bool enable_z,
    typename Real33, typename CoefficientMatrix>
size_t LatticeOperatorInternal(
        const Real33 & abc,
        const size_t ws,
        CoefficientMatrix & lattice,
        const size_t p,
        const size_t maxlayers,
        ssize_t check_convergence_at,
        const double maxerror)
{
    CoefficientMatrix L_star(p);
    L_star_operator<enable_x, enable_y, enable_z>(abc, ws, L_star, p);
    L_star.populate_lower();

    CoefficientMatrix O_star(p);
    O_star_operator<enable_x, enable_y, enable_z>(abc, ws, O_star, p);
    O_star.populate_lower();

    size_t eff_layers = 0;

    if (maxlayers == 0) {
        CopyTriangle(L_star, lattice, p);
    } else {
        CoefficientMatrix L_i(p);
        eff_layers = LatticeOperatorIterate(O_star, L_star, ws, L_i, p,
                maxlayers, check_convergence_at, maxerror);
        CopyTriangle(L_i, lattice, p);
    }
    return eff_layers;
}

// compute lattice operator
// 'lattice' will be overwritten, not accumulated into
template <bool enable_x, bool enable_y, bool enable_z,
         typename Real33, typename CoefficientMatrix>
size_t __attribute__((noinline)) LatticeOperator1D2D(
        const Real33 & abc,
        const size_t ws,
        CoefficientMatrix & lattice,
        const size_t p,
        const size_t maxlayers,
        const double maxerror = std::numeric_limits<double>::min())
{
    //static_assert(enable_x || enable_y || enable_z, "0D-periodic");
    //static_assert(!(enable_x && enable_y && enable_z), "3D-periodic");
    // TODO: generalization to 3D will require computation of the check_convergence_at point
    // and the non-converging elements that are to be zeroed
    size_t layers = LatticeOperatorInternal<enable_x, enable_y, enable_z>(
            abc, ws, lattice, 2 * p, maxlayers, 2, maxerror);
    lattice.zero(1);  // eliminate anything before (2, 0)
    return layers;
}

// compute lattice operator
// 'lattice' will be overwritten, not accumulated into
template <typename Real33, typename CoefficientMatrix>
size_t __attribute__((noinline)) LatticeOperator1D(
        const Real33 & abc,
        const size_t ws,
        CoefficientMatrix & lattice,
        const size_t p,
        const size_t maxlayers,
        const double maxerror = std::numeric_limits<double>::min())
{
    size_t layers = LatticeOperatorInternal<false, false, true>(
            abc, ws, lattice, 2 * p, maxlayers, 2, maxerror);
    lattice.zero(1);  // eliminate anything before (2, 0)
    return layers;
}

// compute lattice operator
// 'lattice' will be overwritten, not accumulated into
template <typename Real33, typename CoefficientMatrix>
size_t __attribute__((noinline)) LatticeOperator2D(
        const Real33 & abc,
        const size_t ws,
        CoefficientMatrix & lattice,
        const size_t p,
        const size_t maxlayers,
        const double maxerror = std::numeric_limits<double>::min())
{
    size_t layers = LatticeOperatorInternal<true, true, false>(
            abc, ws, lattice, 2 * p, maxlayers, 2, maxerror);
    lattice.zero(1);  // eliminate anything before (2, 0)
    return layers;
}

// compute lattice operator
// 'lattice' will be overwritten, not accumulated into
template <typename Real33, typename CoefficientMatrix>
size_t __attribute__((noinline)) LatticeOperator3D(const Real33 & abc,
        const size_t ws,
        CoefficientMatrix & lattice,
        const size_t p,
        const size_t maxlayers,
        const double maxerror/* = std::numeric_limits<double>::min()*/)
{

    ssize_t check_for_convergence = (ssize_t)std::min((int)(2*p),4);
    size_t layers = LatticeOperatorInternal<true, true, true>(
            abc, ws, lattice, 2 * p, maxlayers, check_for_convergence, maxerror);
    lattice.zero(std::max((int)(check_for_convergence-1),0));  // eliminate anything before (4, 0)
    return layers;

}

////////////////////////////////////////////////////////////////////////////

template <typename ESCoefficientMatrix>
__attribute__((noinline))
size_t LatticeOperatorESIterate(
        const ESCoefficientMatrix & O_star,
        const ESCoefficientMatrix & L_star,
        const ssize_t ws,
        ESCoefficientMatrix & L_i,
        const size_t p,
        const size_t maxlayers,
        size_t check_convergence_at = 0,
        const double maxerror = std::numeric_limits<double>::min())
{
    typedef typename ESCoefficientMatrix::value_type::value_type::value_type Real;
    Real oldval = L_star.energy(check_convergence_at, 0).real();
    ESCoefficientMatrix L_pred_scaled(p);
    ESCoefficientMatrix L_tmp0(p);
    ESCoefficientMatrix L_tmp(p);
    for (size_t i = 0; i < maxlayers; ++i) {
        if (i == 0) {
            ScaleL_operator(L_star.energy, ws, L_pred_scaled.energy, p);
            ScaleL_operator(L_star.stress.a, ws, L_pred_scaled.stress.a, p);
            ScaleL_operator(L_star.stress.b, ws, L_pred_scaled.stress.b, p);
            ScaleL_operator(L_star.stress.c, ws, L_pred_scaled.stress.c, p);
        } else {
            ScaleL_operator(L_i.energy, ws, L_pred_scaled.energy, p);
            ScaleL_operator(L_i.stress.a, ws, L_pred_scaled.stress.a, p);
            ScaleL_operator(L_i.stress.b, ws, L_pred_scaled.stress.b, p);
            ScaleL_operator(L_i.stress.c, ws, L_pred_scaled.stress.c, p);
        }
        L_pred_scaled.populate_lower();

        // compute L_{i+1}^{E,a,b,c}

        CopyTriangle(L_star.energy, L_i.energy, p);
        single_rawM2L(O_star.energy, L_pred_scaled.energy, L_i.energy, p);

        CopyTriangle(L_star.stress.a, L_i.stress.a, p);
        L_tmp0.stress.a.zero(p);
        single_rawM2L(O_star.stress.a, L_pred_scaled.energy, L_tmp0.stress.a, p);
        L_tmp.stress.a.zero(p);
        single_rawM2L(O_star.energy, L_pred_scaled.stress.a, L_tmp.stress.a, p);

        CopyTriangle(L_star.stress.b, L_i.stress.b, p);
        L_tmp0.stress.b.zero(p);
        single_rawM2L(O_star.stress.b, L_pred_scaled.energy, L_tmp0.stress.b, p);
        L_tmp.stress.b.zero(p);
        single_rawM2L(O_star.energy, L_pred_scaled.stress.b, L_tmp.stress.b, p);

        CopyTriangle(L_star.stress.c, L_i.stress.c, p);
        L_tmp0.stress.c.zero(p);
        single_rawM2L(O_star.stress.c, L_pred_scaled.energy, L_tmp0.stress.c, p);
        L_tmp.stress.c.zero(p);
        single_rawM2L(O_star.energy, L_pred_scaled.stress.c, L_tmp.stress.c, p);

        Real scale(2 * ws + 1);

        for (ssize_t l = 0; l <= ssize_t(p); ++l) {
            for (ssize_t m = 0; m <= l; ++m) {
                L_i.stress.a(l, m) += L_tmp0.stress.a(l, m) + L_tmp.stress.a(l, m) * scale;
                L_i.stress.b(l, m) += L_tmp0.stress.b(l, m) + L_tmp.stress.b(l, m) * scale;
                L_i.stress.c(l, m) += L_tmp0.stress.c(l, m) + L_tmp.stress.c(l, m) * scale;
            }
        }

        Real newval = L_i.energy(check_convergence_at, 0).real();
        if ((std::abs(oldval - newval) < maxerror))
            return i + 1;
        oldval = newval;
    }
    return maxlayers + 1;
}

template unsigned long LatticeOperator3D<ABC<XYZ<double> >, MultipoleCoefficientsUpper<double, Device<double> > >(ABC<XYZ<double> > const&, unsigned long, MultipoleCoefficientsUpper<double, Device<double> > &, unsigned long, unsigned long, double);
template unsigned long LatticeOperator3D<ABC<XYZ<float> >, MultipoleCoefficientsUpper<float, Device<float> > >(ABC<XYZ<float> > const&, unsigned long, MultipoleCoefficientsUpper<float, Device<float> > &, unsigned long, unsigned long, double);

}//namespace end

#endif
