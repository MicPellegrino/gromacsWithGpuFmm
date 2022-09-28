#ifndef _BK_FMM_particle2local_hpp_
#define _BK_FMM_particle2local_hpp_

#include "global_functions.hpp"
#include "Vec.h"
#include "xyzq.hpp"

namespace gmx_gpu_fmm{

// add local coefficients for a particle with charge q at relative position (x,y,z)
// to the (non-negative triangle) local coefficient matrix upto order p
template <typename Real4, typename CoefficientMatrix>
inline void P2L_old(
        const Real4 & xyzq,
        CoefficientMatrix & mu,
        const size_t p)
{
    printf("running P2L_old?\n");
    typedef typename Real4::value_type real_type;
    typedef typename CoefficientMatrix::value_type complex_type;
    Vec_traits<real_type> VT;

    const real_type x = xyzq.x;
    const real_type y = xyzq.y;
    const real_type z = xyzq.z;
    const real_type q = xyzq.q;

    real_type dist_squared = x * x + y * y + z * z;
    real_type recip_dist = rsqrt(dist_squared);
    real_type recip_dist_squared = reciprocal(dist_squared);

    // mu_0_0  (for the first iteration)
    complex_type mu_mplus1_mplus1(q * recip_dist);

    // mu_0_0 upto mu_p_p-1
    for (size_t m = 0; m < p; ++m) {
        // mu_m_m  (from previous iteration)
        complex_type mu_m_m = mu_mplus1_mplus1;
        mu(m, m) += mu_m_m;

        // mu_m+1_m+1  (for the next iteration)
        mu_mplus1_mplus1 = VT.same(2 * m + 1) * recip_dist_squared * complex_type(x, y) * mu_m_m;

        // mu_m+1_m
        complex_type mu_mplus1_m = VT.same(2 * m + 1) * recip_dist_squared * z * mu_m_m;
        mu(m + 1, m) += mu_mplus1_m;

        // mu_m+2_m upto mu_p_m
        complex_type mu_lminus2_m = mu_m_m;
        complex_type mu_lminus1_m = mu_mplus1_m;
        for (size_t l = m + 2; l <= p; ++l) {
            // mu_l_m
            complex_type mu_l_m = recip_dist_squared *
                (VT.same(2 * l - 1) * z * mu_lminus1_m - VT.same((l - 1) * (l - 1) - m * m) * mu_lminus2_m);
            mu(l, m) += mu_l_m;
            mu_lminus2_m = mu_lminus1_m;
            mu_lminus1_m = mu_l_m;
        }
    }
    // mu_p_p
    mu(p, p) += mu_mplus1_mplus1;
}

// add local coefficients for a particle with charge q at relative position (x,y,z)
// to the (non-negative triangle) local coefficient matrix upto order p
template <typename Real4, typename CoefficientMatrix>
inline void P2L_delay_complex(
        const Real4 & xyzq,
        CoefficientMatrix & mu,
        const size_t p)
{
    printf("running P2L_delay_complex?\n");
    typedef typename Real4::value_type real_type;
    typedef typename CoefficientMatrix::value_type complex_type;
    Vec_traits<real_type> VT;

    const real_type x = xyzq.x;
    const real_type y = xyzq.y;
    const real_type z = xyzq.z;
    const real_type q = xyzq.q;

    real_type dist_squared = x * x + y * y + z * z;
    real_type recip_dist = rsqrt(dist_squared);
    real_type recip_dist_squared = reciprocal(dist_squared);
    real_type z_recip_dist_squared = z * recip_dist_squared;

    complex_type complex_x_y(x, y);

    // mu_0_0  (for the first iteration)
    complex_type complex_part(VT.same(1));
    real_type mu_mplus1_mplus1 = q * recip_dist;

    // mu_0_0 upto mu_p_p-1
    for (size_t m = 0; m < p; ++m) {
        // mu_m_m  (from previous iteration)
        real_type mu_m_m = mu_mplus1_mplus1;
        mu(m, m) += complex_part * mu_m_m;

        // mu_m+1_m+1  (for the next iteration)
        mu_mplus1_mplus1 = VT.same(2 * m + 1) * recip_dist_squared * mu_m_m;

        // mu_m+1_m
        real_type mu_mplus1_m = (VT.same(2 * m + 1) * recip_dist_squared * mu_m_m) * z;
        mu(m + 1, m) += complex_part * mu_mplus1_m;

        // mu_m+2_m upto mu_p_m
        real_type mu_lminus2_m = mu_m_m;
        real_type mu_lminus1_m = mu_mplus1_m;
        for (size_t l = m + 2; l <= p; ++l) {
            // mu_l_m
            real_type mu_l_m =
                VT.same(2 * l - 1) * z_recip_dist_squared * mu_lminus1_m
                - VT.same((l - 1) * (l - 1) - m * m) * recip_dist_squared * mu_lminus2_m;
            mu(l, m) += complex_part * mu_l_m;
            mu_lminus2_m = mu_lminus1_m;
            mu_lminus1_m = mu_l_m;
        }

        // (for the next iteration)
        complex_part *= complex_x_y;
    }
    // mu_p_p
    mu(p, p) += complex_part * mu_mplus1_mplus1;
}


// compute local coefficients upto order p for a particle with charge q
// at relative position (x,y,z) and store (accumulate/assign/...) them in
// a (non-negative triangle) local coefficient matrix
template <typename Real4, typename CoefficientStorage>
inline
__attribute__((flatten))
//__attribute__((noinline))
void genericP2L_reference(
        const Real4 & xyzq,
        CoefficientStorage & mu,
        const size_t p)
{
    printf("running genericP2L_reference?\n");
    typedef typename CoefficientStorage::value_type complex_type;
    typedef typename complex_type::value_type real_type;
    Vec_traits<real_type> VT;

    const real_type x = xyzq.x;
    const real_type y = xyzq.y;
    const real_type z = xyzq.z;
    const real_type q = xyzq.q;

    real_type dist_squared = x * x + y * y + z * z;
    real_type recip_dist = rsqrt(dist_squared);
    real_type recip_dist_squared = reciprocal(dist_squared);
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
        mu(m, m, complex_x_y_m * mu_m_m);

        // mu_m+1_m+1  (for the next iteration)
        mu_mplus1_mplus1 = e_m * mu_m_m;

        // mu_m+1_m
        real_type mu_mplus1_m = mu_mplus1_mplus1 * z;
        mu(m + 1, m, complex_x_y_m * mu_mplus1_m);

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
            mu(l, m, complex_x_y_m * mu_l_m);

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
    mu(p, p, complex_x_y_m * mu_mplus1_mplus1);
}

template <typename T>
T set_min_max(T val)
{
    T result = val;
    /*
    if( val.real() > std::numeric_limits<Real>::max() )
    {
        result.real( std::numeric_limits<Real>::max() );
    }
    if( val.imag() > std::numeric_limits<Real>::max() )
    {
        result.imag( std::numeric_limits<Real>::max() );
    }
    if( val.real() < std::numeric_limits<Real>::min() )
    {
        result.real( std::numeric_limits<Real>::min() );
    }
    if( val.imag() < std::numeric_limits<Real>::min() )
    {
        result.imag( std::numeric_limits<Real>::min() );
    }
    */
    return result;
}


// add local coefficients for a particle with charge q at relative position (x,y,z)
// to the (non-negative triangle) local coefficient matrix upto order p
template <typename Real4, typename CoefficientMatrix>
inline
__attribute__((flatten))
void P2L(
        const Real4 & xyzq,
        CoefficientMatrix & mu,
        const size_t p)
{
    //printf("running p2l?\n");
    typedef typename Real4::value_type real_type;
    typedef typename CoefficientMatrix::value_type complex_type;

    //typedef double real_type;
    //typedef fmm::complex<real_type> complex_type;
    Vec_traits<real_type> VT;

    const real_type x = xyzq.x;
    const real_type y = xyzq.y;
    const real_type z = xyzq.z;
    const real_type q = xyzq.q;

    real_type dist_squared = x * x + y * y + z * z;
    real_type recip_dist = rsqrt(dist_squared);
    real_type recip_dist_squared = reciprocal(dist_squared);
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
        mu(m, m) += set_min_max(complex_x_y_m * mu_m_m);
        //printf("%e %e\n", set_min_max(complex_x_y_m * mu_m_m).real(), set_min_max(complex_x_y_m * mu_m_m).imag() );

        // mu_m+1_m+1  (for the next iteration)
        mu_mplus1_mplus1 = e_m * mu_m_m;

        // mu_m+1_m
        real_type mu_mplus1_m = mu_mplus1_mplus1 * z;
        mu(m + 1, m) += set_min_max(complex_x_y_m * mu_mplus1_m);
        //printf("%e %e\n",(complex_x_y_m * mu_mplus1_m).real(),(complex_x_y_m * mu_mplus1_m).imag());

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
            mu(l, m) += set_min_max(complex_x_y_m * mu_l_m);
            //printf("%e %e\n",(complex_x_y_m * mu_l_m).real(),(complex_x_y_m * mu_l_m).imag());

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
    //printf("%e %e\n",(complex_x_y_m * mu_mplus1_mplus1).real(),(complex_x_y_m * mu_mplus1_mplus1).imag());
    mu(p, p) += set_min_max(complex_x_y_m * mu_mplus1_mplus1);
}

// add local coefficients for a particle with charge q at relative position (x,y,z)
// to the (non-negative triangle) local coefficient matrix upto order p
template <typename Real4, typename CoefficientStorage>
inline
//__attribute__((noinline))
__attribute__((flatten))
void genericP2L_rowit(
        const Real4 & xyzq,
        CoefficientStorage & mu,
        const size_t p)
{
    typedef typename CoefficientStorage::row_iterator row_iterator;
    typedef typename CoefficientStorage::value_type complex_type;
    typedef typename complex_type::value_type real_type;
    Vec_traits<real_type> VT;

    const real_type x = xyzq.x;
    const real_type y = xyzq.y;
    const real_type z = xyzq.z;
    const real_type q = xyzq.q;

    real_type dist_squared = x * x + y * y + z * z;
    real_type recip_dist = rsqrt(dist_squared);
    real_type recip_dist_squared = reciprocal(dist_squared);
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
        row_iterator mu_row_m = mu.get_row_iterator(m, m);
        mu_row_m.prefetch(0 * 64);
        mu_row_m.prefetch(4 * 64);

        // mu_m_m  (from previous iteration)
        real_type mu_m_m = mu_mplus1_mplus1;
        (*mu_row_m)(complex_x_y_m * mu_m_m);
        ++mu_row_m;

        // mu_m+1_m+1  (for the next iteration)
        mu_mplus1_mplus1 = e_m * mu_m_m;

        // mu_m+1_m
        real_type mu_mplus1_m = mu_mplus1_mplus1 * z;
        (*mu_row_m)(complex_x_y_m * mu_mplus1_m);
        ++mu_row_m;

        // mu_m+2_m upto mu_p_m
        real_type e_mplus1 = e_m + twice_recip_dist_squared;
        real_type mu_lminus2_m = mu_m_m;
        real_type mu_lminus1_m = mu_mplus1_m;
        real_type f_l = e_mplus1 * z;           // iterative computation of ((2*l - 1)*z/R^2)
        real_type h_l = e_m;                    // iterative computation of ((2*l - 3)/R^2)
        real_type g_l = h_l;                    // iterative computation of (((l - 1)^2 - m^2)/R^2)
        for (size_t l = m + 2; l <= p; ++l) {
            mu_row_m.prefetch(8 * 64);

            // mu_l_m
            real_type mu_l_m = f_l * mu_lminus1_m - g_l * mu_lminus2_m;
            (*mu_row_m)(complex_x_y_m * mu_l_m);
            ++mu_row_m;

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
    row_iterator mu_p_p = mu.get_row_iterator(p, p);
    (*mu_p_p)(complex_x_y_m * mu_mplus1_mplus1);
}


// DEPRECATED name
template <typename Real4, typename CoefficientStorage>
inline
void P2L_generic(
        const Real4 & xyzq,
        CoefficientStorage & mu,
        const size_t p)
{
    genericP2L_reference(xyzq, mu, p);
}


// Particle2Local
// add coefficients contributed by the given particles to mu
// (input mu must be either zeroed or contain proper coefficients to be added to)
template <typename TableReal, typename Real3, typename CoefficientMatrix>
void P2L(
        const TableReal * vx,
        const TableReal * vy,
        const TableReal * vz,
        const TableReal * vq,
        const Real3 & expansion_point,
        const size_t begin,
        const size_t end,
        CoefficientMatrix & mu,
        const size_t p)
{
    typedef typename Real3::value_type Real;

    for (size_t particle = begin; particle < end; ++particle) {
        Real3 xyz(vx[particle], vy[particle], vz[particle]);
        Real q = vq[particle];
        P2L(make_xyzq(xyz - expansion_point, q), mu, p);
    }
}

}//namespace end

#endif
