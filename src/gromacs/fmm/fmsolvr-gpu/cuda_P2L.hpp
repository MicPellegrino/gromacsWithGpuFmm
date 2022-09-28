#ifndef _BK_FMM_cuda_P2L_hpp
#define _BK_FMM_cuda_P2L_hpp

#include "cub_h.h"
#include "cuda_keywords.hpp"

namespace gmx_gpu_fmm{

namespace p2l{
template <typename Real>
DEVICE
inline Real __reciprocal(Real a)
{
    return Real(1.) / a;
}

}

// add local coefficients for a particle with charge q at relative position (x,y,z)
// to the (non-negative triangle) local coefficient matrix upto order p
template <typename Real4, typename CoeffMatrix>
__device__
void __P2L(
        const Real4 & xyzq,
        CoeffMatrix & mu,
        const size_t p)
{
    typedef typename Real4::value_type real_type;
    typedef typename CoeffMatrix::value_type complex_type;
    //typedef double real_type;
    //typedef fmm::complex<real_type> complex_type;
    Vec_traits<real_type> VT;

    const real_type x = xyzq.x;
    const real_type y = xyzq.y;
    const real_type z = xyzq.z;
    const real_type q = xyzq.q;

    real_type dist_squared = x * x + y * y + z * z;
    real_type recip_dist = ::rsqrt(dist_squared);
    real_type recip_dist_squared = p2l::__reciprocal(dist_squared);
    real_type twice_recip_dist_squared = recip_dist_squared + recip_dist_squared;
    real_type z_recip_dist_squared = z * recip_dist_squared;
    real_type twice_z_recip_dist_squared = z_recip_dist_squared + z_recip_dist_squared;

    complex_type complex_x_y(x, y);
    complex_type complex_x_y_m(VT.same(1));     // iterative computation of ((x + iy)^m)

    // mu_0_0  (for the first iteration)
    real_type mu_mplus1_mplus1 = q * recip_dist;
    real_type e_m = recip_dist_squared;         // iterative computation of ((2*m + 1)/R^2)

    // mu_0_0 upto mu_p_p-1
    for (size_t m = 0; m < p; ++m) {
        // mu_m_m  (from previous iteration)
        real_type mu_m_m = mu_mplus1_mplus1;
        mu(m, m) %= complex_x_y_m * mu_m_m;

        // mu_m+1_m+1  (for the next iteration)
        mu_mplus1_mplus1 = e_m * mu_m_m;

        // mu_m+1_m
        real_type mu_mplus1_m = mu_mplus1_mplus1 * z;
        mu(m + 1, m) %= complex_x_y_m * mu_mplus1_m;

        // mu_m+2_m upto mu_p_m
        real_type e_mplus1 = e_m + twice_recip_dist_squared;
        real_type mu_lminus2_m = mu_m_m;
        real_type mu_lminus1_m = mu_mplus1_m;
        real_type f_l = e_mplus1 * z;           // iterative computation of ((2*l - 1)*z/R^2)
        real_type h_l = e_m;                    // iterative computation of ((2*l - 3)/R^2)
        real_type g_l = h_l;                    // iterative computation of (((l - 1)^2 - m^2)/R^2)
        for (size_t l = m + 2; l <= p; ++l) {
            // mu_l_m
            real_type mu_l_m = f_l * mu_lminus1_m - g_l * mu_lminus2_m;
            mu(l, m) %= complex_x_y_m * mu_l_m;

            // (for the next iteration)
            mu_lminus2_m = mu_lminus1_m;
            mu_lminus1_m = mu_l_m;
            f_l += twice_z_recip_dist_squared;
            h_l += twice_recip_dist_squared;
            g_l += h_l;
        }

        // (for the next iteration)
        complex_x_y_m *= complex_x_y;
        e_m = e_mplus1;
    }

    // mu_p_p
    mu(p, p) %= complex_x_y_m * mu_mplus1_mplus1;
}

// add local coefficients for a particle with charge q at relative position (x,y,z)
// to the (non-negative triangle) local coefficient matrix upto order p
template <typename Real4, typename CoeffMatrix>
__device__
void __P2L_cubereduce(
        const Real4 & xyzq,
        CoeffMatrix & mu,
        const size_t p)
{

    typedef typename Real4::value_type real_type;
    typedef typename CoeffMatrix::value_type complex_type;
    //typedef double real_type;
    //typedef complex<real_type> complex_type;
    Vec_traits<real_type> VT;

    typedef cub::WarpReduce<complex_type> WarpReducer;
    __shared__ typename WarpReducer::TempStorage temp[128];

    const real_type x = xyzq.x;
    const real_type y = xyzq.y;
    const real_type z = xyzq.z;
    const real_type q = xyzq.q;

    real_type dist_squared = x * x + y * y + z * z;
    real_type recip_dist = ::rsqrt(dist_squared);
    real_type recip_dist_squared = p2l::__reciprocal(dist_squared);
    real_type twice_recip_dist_squared = recip_dist_squared + recip_dist_squared;
    real_type z_recip_dist_squared = z * recip_dist_squared;
    real_type twice_z_recip_dist_squared = z_recip_dist_squared + z_recip_dist_squared;

    complex_type complex_x_y(x, y);
    complex_type complex_x_y_m(VT.same(1));     // iterative computation of ((x + iy)^m)

    // mu_0_0  (for the first iteration)
    real_type mu_mplus1_mplus1 = q * recip_dist;
    real_type e_m = recip_dist_squared;         // iterative computation of ((2*m + 1)/R^2)

    size_t lane_id = threadIdx.x%32;
    // mu_0_0 upto mu_p_p-1
    for (size_t m = 0; m < p; ++m) {
        // mu_m_m  (from previous iteration)
        real_type mu_m_m = mu_mplus1_mplus1;

        complex_type out_mu_m_m = complex_x_y_m * mu_m_m;
        complex_type sum_out_mu_m_m = WarpReducer(temp[lane_id]).Sum(out_mu_m_m);
        if(lane_id == 0)
            mu(m, m) %=sum_out_mu_m_m;

        // mu_m+1_m+1  (for the next iteration)
        mu_mplus1_mplus1 = e_m * mu_m_m;

        // mu_m+1_m
        real_type mu_mplus1_m = mu_mplus1_mplus1 * z;

        complex_type out_mu_mplus1_m = complex_x_y_m * mu_mplus1_m;
        complex_type sum_out_mu_mplus1_m = WarpReducer(temp[lane_id]).Sum(out_mu_mplus1_m);
        if(lane_id == 0)
            mu(m + 1, m) %= sum_out_mu_mplus1_m;

        // mu_m+2_m upto mu_p_m
        real_type e_mplus1 = e_m + twice_recip_dist_squared;
        real_type mu_lminus2_m = mu_m_m;
        real_type mu_lminus1_m = mu_mplus1_m;
        real_type f_l = e_mplus1 * z;           // iterative computation of ((2*l - 1)*z/R^2)
        real_type h_l = e_m;                    // iterative computation of ((2*l - 3)/R^2)
        real_type g_l = h_l;                    // iterative computation of (((l - 1)^2 - m^2)/R^2)
        for (size_t l = m + 2; l <= p; ++l) {
            // mu_l_m
            real_type mu_l_m = f_l * mu_lminus1_m - g_l * mu_lminus2_m;

            complex_type out_mu_l_m = complex_x_y_m * mu_l_m;
            complex_type sum_out_mu_l_m = WarpReducer(temp[lane_id]).Sum(out_mu_l_m);
            if(lane_id == 0)
                mu(l, m) %= sum_out_mu_l_m;

            // (for the next iteration)
            mu_lminus2_m = mu_lminus1_m;
            mu_lminus1_m = mu_l_m;
            f_l += twice_z_recip_dist_squared;
            h_l += twice_recip_dist_squared;
            g_l += h_l;
        }

        // (for the next iteration)
        complex_x_y_m *= complex_x_y;
        e_m = e_mplus1;
    }

    // mu_p_p
    complex_type out_mu_p_p = complex_x_y_m * mu_mplus1_mplus1;
    complex_type sum_out_mu_p_p = WarpReducer(temp[lane_id]).Sum(out_mu_p_p);
    if(lane_id == 0)
        mu(p, p) %=sum_out_mu_p_p;
}

}//namespace end

#endif
