#ifndef _BK_FMM_lstar_operator_hpp_
#define _BK_FMM_lstar_operator_hpp_

#include <iomanip>

#include "particle2local.hpp"
#include "Vec.h"
#include "triangular_array.hpp"
#include "generic_accumulate.hpp"
#include "super_accumulate.hpp"
#include "iterables.hpp"
#include "triangular_array_adaptor.hpp"
#include "energy_stress_triangular_array_group.hpp"

// for benchmarking and disassembling
#define lstar_layersum_noinline __attribute__((noinline))




namespace gmx_gpu_fmm{

//
// L_* computation is an *assigning* operation.
// The result will be placed in lstar (a reference to some triangular
// array container), either by copying it there from some temporary
// storage or by first zeroing it and thereafter accumulating into it.
// lattice contribution from the box centered at (x,y,z)
template <bool compute_complex, typename Real, typename CoefficientMatrix>
inline void Lattice3Dcubic(
        const Real x,
        const Real y,
        const Real z,
        const Real s,
        CoefficientMatrix & lattice,
        const size_t p_)
{
    typedef Real real_type;
    typedef typename CoefficientMatrix::complex_type complex_type;

    const size_t p = p_ - (p_ & 1);  // ensure p is even
    real_type dist_squared = x * x + y * y + z * z;
    real_type recip_dist = rsqrt(dist_squared);
    real_type recip_dist_squared = reciprocal(dist_squared);
    real_type recip_dist_to_the_4 = recip_dist_squared * recip_dist_squared;
    complex_type x_y_squared;
    if (compute_complex)
        x_y_squared = complex_type(x * x - y * y, x * y + x * y);  // (x + iy)^2
    else
        x_y_squared = complex_type(x * x - y * y);

    // lattice_0_0  (for the first iteration)
    complex_type lattice_mplus2_mplus2(recip_dist * s);

    // lattice_0_0 upto lattice_p_p-1
    for (size_t m = 0; m < p; m += 2) {
        // lattice_m_m  (from previous iteration)
        complex_type lattice_m_m = lattice_mplus2_mplus2;
        lattice(m, m) += lattice_m_m;

        // lattice_m+2_m+2  (for the next iteration)
        lattice_mplus2_mplus2 = real_type((2 * m + 3) * (2 * m + 1)) * recip_dist_to_the_4 * x_y_squared * lattice_m_m;
        //if (m == 4)
        //    std::cout << x << " " << y << " " << z << " " << lattice_mplus2_mplus2 << std::endl;

        // lattice_m+1_m
        complex_type lattice_mplus1_m = real_type(2 * m + 1) * recip_dist_squared * z * lattice_m_m;
        // don't store

        // lattice_m+2_m upto lattice_p_m
        complex_type lattice_lminus2_m = lattice_m_m;
        complex_type lattice_lminus1_m = lattice_mplus1_m;
        for (size_t l = m + 2; l <= p; l += 2) {
            // lattice_l_m
            complex_type lattice_l_m = recip_dist_squared *
                (real_type(2 * l - 1) * z * lattice_lminus1_m - real_type((l - 1) * (l - 1) - m * m) * lattice_lminus2_m);
            lattice(l, m) += lattice_l_m;

            complex_type lattice_lplus1_m = recip_dist_squared *
                (real_type(2 * l + 1) * z * lattice_l_m - real_type((l) * (l) - m * m) * lattice_lminus1_m);
            // don't store

            lattice_lminus2_m = lattice_l_m;
            lattice_lminus1_m = lattice_lplus1_m;

            //if (l == 4 && m == 0)
            //    std::cout << x << " " << y << " " << z << " " << lattice_l_m << std::endl;
        }
    }

    // lattice_p_p
    lattice(p, p) += lattice_mplus2_mplus2;
}

// sum up (-outer, -outer, -outer) .. (outer, outer, outer)
// excluding (-inner, -inner, -inner) .. (inner, inner, inner)
template <bool enable_x, bool enable_y, bool enable_z,
    typename Real33, typename CoefficientMatrix>
size_t Lattice3Dcubic_layer_sum(
        const Real33 & abc,
        const ssize_t inner,
        const ssize_t outer,
        CoefficientMatrix & accu,
        const size_t p)
{
    typedef typename Real33::value_type Real3;
    typedef typename Real3::value_type Real;

    const bool do_complex = true;
    const ssize_t outer_x = outer;
    const ssize_t outer_y = outer;
    const ssize_t outer_z = outer;

    size_t cnt = 0;
    for (ssize_t i = do_complex ? -outer_x : 0; i <= +outer_x; ++i)
        for (ssize_t j = do_complex ? -outer_z : 0; j <= +outer_y; ++j)
            for (ssize_t k = do_complex ? -outer_z : 0; k <= +outer_z; ++k)
                if (i < -inner || i > +inner || j < -inner || j > +inner || k < -inner || k > +inner) {
                    ++cnt;
                    Real3 m = abc.a * i + abc.b * j + abc.c * k;
                    if (!do_complex) {
                        size_t zeroes = (i == 0) + (j == 0) + (k == 0);
                        //Real s = zeroes == 0 ? 1.0 : zeroes == 1 ? 0.5 : zeroes == 2 ? 0.25 : -0.;
                        Real s = zeroes == 0 ? 8.0 : zeroes == 1 ? 4.0 : zeroes == 2 ? 2.0 : -0.;
                        //std::cout << i << " " << j << " " << k << " " << s << std::endl;
                        Lattice3Dcubic<false>(m.x, m.y, m.z, s, accu, p);
                    } else {
                        Lattice3Dcubic<true>(m.x, m.y, m.z, Real(1.0), accu, p);
                    }
                }
    return cnt;
}


// sum up (-outer, -outer, -outer) .. (outer, outer, outer)
// excluding (-inner, -inner, -inner) .. (inner, inner, inner)
template <bool enable_a, bool enable_b, bool enable_c,
    typename Real33, typename CoefficientMatrix>
lstar_layersum_noinline
size_t L_star_layer_sum_reference(
        const Real33 & abc,
        const ssize_t inner,
        const ssize_t outer,
        CoefficientMatrix & lstar,
        const size_t p)
{
    typedef typename Real33::value_type Real3;
    typedef typename Real3::value_type Real;

    assert(p <= lstar.p());
    lstar.zero(p);

    const ssize_t outer_a = enable_a ? outer : 0;
    const ssize_t outer_b = enable_b ? outer : 0;
    const ssize_t outer_c = enable_c ? outer : 0;

    size_t cnt = 0;
    for (ssize_t i = -outer_a; i <= +outer_a; ++i)
        for (ssize_t j = -outer_b; j <= +outer_b; ++j)
            for (ssize_t k = -outer_c; k <= +outer_c; ++k)
                if (i < -inner || i > +inner || j < -inner || j > +inner || k < -inner || k > +inner) {
                    ++cnt;
                    Real3 m = abc * Real3(i, j, k);
                    //std::cout<<m<<std::endl;
                    P2L(make_xyzq(m, Real(1.0)), lstar, p);  // accumulate
                }
    return cnt;
}


template <bool enable_x, bool enable_y, bool enable_z,
    typename Real33, typename CoefficientMatrix>
size_t L_star_cube_sum_reference_dump(
        const Real33 & abc,
        const ssize_t inner,
        const ssize_t outer,
        CoefficientMatrix & lstar,
        const size_t p)
{
    typedef typename Real33::value_type Real3;
    typedef typename Real3::value_type Real;

    assert(p <= lstar.p());
    lstar.zero(p);

    const ssize_t outer_x = enable_x ? outer : 0;
    const ssize_t outer_y = enable_y ? outer : 0;
    const ssize_t outer_z = enable_z ? outer : 0;

    std::cout << std::scientific << std::setprecision(30);
    std::cout << "# i\tj\tk\tl\tm\treal\timag" << std::endl;

    size_t cnt = 0;
    for (ssize_t i = -outer_x; i <= +outer_x; ++i)
        for (ssize_t j = -outer_y; j <= +outer_y; ++j)
            for (ssize_t k = -outer_z; k <= +outer_z; ++k)
                if (i < -inner || i > +inner || j < -inner || j > +inner || k < -inner || k > +inner) {
                    ++cnt;
                    Real3 m = abc.a * i + abc.b * j + abc.c * k;
                    lstar.zero();
                    P2L(make_xyzq(m.x, m.y, m.z, Real(1.0)), lstar, p);
                    for (ssize_t l = 0; l <= ssize_t(p); ++l)
                        for (ssize_t m = 0; m <= l; ++m)
                            std::cout << i << "\t" << j << "\t" << k << "\t" << l << "\t" << m << "\t" << lstar(l, m).real() << "\t" << lstar(l, m).imag() << std::endl;
                }
    return cnt;
}


// sum up (-outer, -outer, -outer) .. (outer, outer, outer)
// excluding (-inner, -inner, -inner) .. (inner, inner, inner)
template <bool enable_a, bool enable_b, bool enable_c,
    typename Real33, typename CoefficientMatrix>
lstar_layersum_noinline
size_t L_star_layer_sum_generic_colstore(
        const Real33 & abc,
        const ssize_t inner,
        const ssize_t outer,
        CoefficientMatrix & lstar,
        const size_t p)
{
    typedef typename Real33::value_type Real3;
    typedef typename Real3::value_type Real;
    typedef typename CoefficientMatrix::value_type complex_type;

    typedef simple_triangular_array_adaptor<CoefficientMatrix, complex_type, just_accumulate> storage_adapter_type;

    assert(p <= lstar.p());
    lstar.zero(p);

    const ssize_t outer_a = enable_a ? outer : 0;
    const ssize_t outer_b = enable_b ? outer : 0;
    const ssize_t outer_c = enable_c ? outer : 0;

    storage_adapter_type storage_adapter(lstar);

    size_t cnt = 0;
    for (ssize_t i = -outer_a; i <= +outer_a; ++i) {
        for (ssize_t j = -outer_b; j <= +outer_b; ++j) {
            for (ssize_t k = -outer_c; k <= +outer_c; ++k) {
                if (i < -inner || i > +inner || j < -inner || j > +inner || k < -inner || k > +inner) {
                    ++cnt;
                    Real3 m = abc * Real3(i, j, k);
                    P2L_generic(make_xyzq(m, Real(1.0)), storage_adapter, p);
                }
            }
        }
    }
    return cnt;
}


// sum up (-outer, -outer, -outer) .. (outer, outer, outer)
// excluding (-inner, -inner, -inner) .. (inner, inner, inner)
template <bool enable_a, bool enable_b, bool enable_c,
    typename Real33, typename CoefficientMatrix>
lstar_layersum_noinline
size_t L_star_layer_sum_generic_rowstore(
        const Real33 & abc,
        const ssize_t inner,
        const ssize_t outer,
        CoefficientMatrix & lstar,
        const size_t p)
{
    typedef typename Real33::value_type Real3;
    typedef typename Real3::value_type Real;
    typedef typename CoefficientMatrix::value_type complex_type;

    typedef UpperTriangularArrayRowMajor<complex_type, void> temp_triangle_type;
    typedef simple_triangular_array_adaptor<temp_triangle_type, complex_type, just_accumulate> storage_adapter_type;

    assert(p <= lstar.p());

    const ssize_t outer_a = enable_a ? outer : 0;
    const ssize_t outer_b = enable_b ? outer : 0;
    const ssize_t outer_c = enable_c ? outer : 0;

    temp_triangle_type temp_storage(p);
    storage_adapter_type storage_adapter(temp_storage);

    size_t cnt = 0;
    for (ssize_t i = -outer_a; i <= +outer_a; ++i) {
        for (ssize_t j = -outer_b; j <= +outer_b; ++j) {
            for (ssize_t k = -outer_c; k <= +outer_c; ++k) {
                if (i < -inner || i > +inner || j < -inner || j > +inner || k < -inner || k > +inner) {
                    ++cnt;
                    Real3 m = abc * Real3(i, j, k);
                    P2L_generic(make_xyzq(m, Real(1.0)), storage_adapter, p);
                }
            }
        }
    }
    for (ssize_t l = 0; l <= ssize_t(p); ++l)
        for (ssize_t m = 0; m <= l; ++m)
            lstar(l, m) = temp_storage(l, m);
    return cnt;
}


// sum up (-outer, -outer, -outer) .. (outer, outer, outer)
// excluding (-inner, -inner, -inner) .. (inner, inner, inner)
template <bool enable_a, bool enable_b, bool enable_c,
    typename Real33, typename ESCoefficientMatrix>
lstar_layersum_noinline
size_t L_star_es_layer_sum_generic_colstore(
        const Real33 & abc,
        const ssize_t inner,
        const ssize_t outer,
        ESCoefficientMatrix & lstar_es,
        const size_t p)
{
    typedef typename Real33::value_type Real3;
    typedef typename Real3::value_type Real;
    typedef typename ESCoefficientMatrix::value_type::value_type complex_type;

    typedef EnergyStressStorageOperator<ESCoefficientMatrix, complex_type, just_accumulate> store_energy_stress_operator_type;
    typedef custom_triangular_array_adaptor<store_energy_stress_operator_type> storage_adapter_type;

    assert(p <= lstar_es.energy.p());
    lstar_es.zero(p);

    const ssize_t outer_a = enable_a ? outer : 0;
    const ssize_t outer_b = enable_b ? outer : 0;
    const ssize_t outer_c = enable_c ? outer : 0;

    size_t cnt = 0;
    for (ssize_t i = -outer_a; i <= +outer_a; ++i) {
        for (ssize_t j = -outer_b; j <= +outer_b; ++j) {
            for (ssize_t k = -outer_c; k <= +outer_c; ++k) {
                if (i < -inner || i > +inner || j < -inner || j > +inner || k < -inner || k > +inner) {
                    ++cnt;
                    Real3 ijk(i, j, k);
                    Real3 m = abc * ijk;
                    storage_adapter_type storage_adapter(lstar_es, ijk.x, ijk.y, ijk.z);
                    P2L_generic(make_xyzq(m, Real(1.0)), storage_adapter, p);
                }
            }
        }
    }
    return cnt;
}


// sum up (-outer, -outer, -outer) .. (outer, outer, outer)
// excluding (-inner, -inner, -inner) .. (inner, inner, inner)
template <bool enable_a, bool enable_b, bool enable_c,
    typename Real33, typename ESCoefficientMatrix>
lstar_layersum_noinline
size_t L_star_es_layer_sum_generic_rowstore(
        const Real33 & abc,
        const ssize_t inner,
        const ssize_t outer,
        ESCoefficientMatrix & lstar_es,
        const size_t p)
{
    typedef typename Real33::value_type Real3;
    typedef typename Real3::value_type Real;
    typedef typename ESCoefficientMatrix::value_type::value_type complex_type;

    typedef UpperTriangularArrayRowMajor<complex_type, void> temp_triangular_type;
    typedef EnergyStressTriangularArrayGroup<temp_triangular_type> es_temp_triangular_type;
    typedef EnergyStressStorageOperator<es_temp_triangular_type, complex_type, just_accumulate> store_energy_stress_operator_type;
    typedef custom_triangular_array_adaptor<store_energy_stress_operator_type> storage_adapter_type;

    assert(p <= lstar_es.energy.p());

    const ssize_t outer_a = enable_a ? outer : 0;
    const ssize_t outer_b = enable_b ? outer : 0;
    const ssize_t outer_c = enable_c ? outer : 0;

    es_temp_triangular_type temp_storage(p);

    size_t cnt = 0;
    for (ssize_t i = -outer_a; i <= +outer_a; ++i) {
        for (ssize_t j = -outer_b; j <= +outer_b; ++j) {
            for (ssize_t k = -outer_c; k <= +outer_c; ++k) {
                if (i < -inner || i > +inner || j < -inner || j > +inner || k < -inner || k > +inner) {
                    ++cnt;
                    Real3 ijk(i, j, k);
                    Real3 m = abc * ijk;
                    storage_adapter_type storage_adapter(temp_storage, ijk.x, ijk.y, ijk.z);
                    P2L_generic(make_xyzq(m, Real(1.0)), storage_adapter, p);
                }
            }
        }
    }
    for (ssize_t l = 0; l <= ssize_t(p); ++l)
        for (ssize_t m = 0; m <= l; ++m) {
            lstar_es.energy(l, m) = temp_storage.energy(l, m);
            lstar_es.stress.a(l, m) = temp_storage.stress.a(l, m);
            lstar_es.stress.b(l, m) = temp_storage.stress.b(l, m);
            lstar_es.stress.c(l, m) = temp_storage.stress.c(l, m);
        }
    return cnt;
}


template <typename Real3>
inline typename Real3::value_type len_squared(const Real3 & xyz)
{
    return xyz.x * xyz.x + xyz.y * xyz.y + xyz.z * xyz.z;
}

template <typename Real3>
inline typename Real3::value_type dot_prod(const Real3 & a, const Real3 & b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

template <typename Real3>
inline Real3 cross_prod(const Real3 & a, const Real3 & b)
{
    return Real3(
            a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x);
}

// abc is a parallelepiped given by three vectors
template <typename Real33>
typename Real33::value_type::value_type umkugeldiameter_squared(const Real33 & abc)
{
    typedef typename Real33::value_type Real3;
    typedef typename Real3::value_type Real;
    Real3 a = abc.a;
    Real3 b = abc.b;
    Real3 c = abc.c;

    // 4 diagonals
    Real l0 = len_squared(a + b + c);
    Real l1 = len_squared(a + b - c);
    Real l2 = len_squared(a + c - b);
    Real l3 = len_squared(b + c - a);

    return std::max(std::max(l0, l1), std::max(l2, l3));
}

// distance of point A from the plane defined by vectors b,c (and origin 0,0,0)
template <typename Real3>
typename Real3::value_type distance_squared(const Real3 & a, const Real3 & b, const Real3 & c)
{
    typedef typename Real3::value_type Real;

    Real3 bxc = cross_prod(b, c);
    Real a_bxc = dot_prod(a, bxc);
    return a_bxc * a_bxc / len_squared(bxc);
}

// abc is a parallelepiped given by three vectors
template <typename Real33>
typename Real33::value_type::value_type inkugeldiameter_squared(const Real33 & abc)
{
    typedef typename Real33::value_type Real3;
    typedef typename Real3::value_type Real;
    Real3 a = abc.a;
    Real3 b = abc.b;
    Real3 c = abc.c;

    // 3 heights
    Real ha = distance_squared(a, b, c);
    Real hb = distance_squared(b, a, c);
    Real hc = distance_squared(c, a, b);

    return std::min(std::min(ha, hb), hc);
}

inline double r2component(ssize_t i)
{
    double d = (i == 0) ? 0. : double(std::abs(i)) - .5;
    return d * d;
}

// sum up (-outer, -outer, -outer) .. (outer, outer, outer)
// excluding (-inner, -inner, -inner) .. (inner, inner, inner)
template <bool enable_x, bool enable_y, bool enable_z,
    typename Real33, typename CoefficientMatrix>
size_t L_star_layer_sum_adaptive(
        const Real33 & abc,
        const ssize_t inner,
        const ssize_t outer,
        CoefficientMatrix & accu,
        const size_t p_max)
{
    typedef typename Real33::value_type Real3;
    typedef typename Real3::value_type Real;

    const ssize_t outer_x = enable_x ? outer : 0;
    const ssize_t outer_y = enable_y ? outer : 0;
    const ssize_t outer_z = enable_z ? outer : 0;

    double a = umkugeldiameter_squared(abc) * .25;
    Real innerscale = 2 * inner + 1;
    double r1 = inkugeldiameter_squared(Real33(abc.a * innerscale, abc.b * innerscale, abc.c * innerscale)) * .25;
    double loga = std::log(a);
    double logr1 = std::log(r1);

    size_t cnt = 0;
    for (ssize_t i = -outer_x; i <= +outer_x; ++i)
        for (ssize_t j = -outer_y; j <= +outer_y; ++j)
            for (ssize_t k = -outer_z; k <= +outer_z; ++k)
                if (i < -inner || i > +inner || j < -inner || j > +inner || k < -inner || k > +inner) {
                    ++cnt;
                    Real3 ijk(i, j, k);
                    Real3 m = abc * ijk;
                    double r22 = r2component(i) + r2component(j) + r2component(k);
                    size_t p = std::ceil(double(p_max) * (loga - logr1) / (loga - std::log(r22)));
                    P2L(make_xyzq(m.x, m.y, m.z, Real(1.0)), accu, p);
                }
    return cnt;
}

#if 0
//#ifdef __AVX__
// sum up (-outer, -outer, -outer) .. (outer, outer, outer)
// excluding (-inner, -inner, -inner) .. (inner, inner, inner)
template <bool enable_x, bool enable_y, bool enable_z,
    typename Real33, typename CoefficientMatrix>
inline size_t L_star_layer_sum_vec(
        const Real33 & abc,
        const ssize_t inner,
        const ssize_t outer,
        CoefficientMatrix & accu,
        const size_t p)
{
    typedef typename Real33::value_type Real3;
    typedef typename Real3::value_type Real;
    typedef typename CoefficientMatrix::value_type complex_type;

    typedef typename simd_engine_default<Real>::type simd_type;
#if 0
    typedef simd_type RealVec;
#else
    typedef SuperSIMD<simd_type, 2> super_simd_type;
    typedef super_simd_type RealVec;
#endif
#define ACCVARIANT 1
#if ACCVARIANT == 0
    // compute with SuperSIMD, temporary storage: SuperSIMD
    typedef RealVec temp_storage_type;
#elif ACCVARIANT == 1
    // compute with SuperSIMD, temporary storage: SIMD
    typedef simd_type temp_storage_type;
#endif
    typedef typename COMPLEX_GENERATOR<RealVec>::type simd_complex_type;
    typedef typename COMPLEX_GENERATOR<temp_storage_type>::type temp_complex_type;
    typedef aligned_allocator<temp_complex_type> temp_allocator_type;
    typedef UpperTriangularArrayRowMajor<temp_complex_type, void, temp_allocator_type> temp_triangle_type;
#if ACCVARIANT == 0
    typedef row_iterable<temp_triangle_type, simd_complex_type, just_accumulate> temp_accumulator_type;
#elif ACCVARIANT == 1
    typedef row_iterable<temp_triangle_type, simd_complex_type, just_complex_super_accumulate> temp_accumulator_type;
#endif

    const ssize_t outer_x = enable_x ? outer : 0;
    const ssize_t outer_y = enable_y ? outer : 0;
    const ssize_t outer_z = enable_z ? outer : 0;

    Vec_traits<RealVec> VT;
    RealVec rvi, rvj, rvk;
    size_t vi = 0;
    temp_triangle_type accu_vec(p);
    temp_accumulator_type temp_accu(accu_vec);

    size_t cnt = 0;
    for (ssize_t i = -outer_x; i <= +outer_x; ++i) {
        Real ri(i);
        for (ssize_t j = -outer_y; j <= +outer_y; ++j) {
            Real rj(j);
            for (ssize_t k = -outer_z; k <= +outer_z; ++k) {
                if (i < -inner || i > +inner || j < -inner || j > +inner || k < -inner || k > +inner) {
                    Real rk(k);
                    ++cnt;
                    VT.element(rvi, vi) = ri;
                    VT.element(rvj, vi) = rj;
                    VT.element(rvk, vi) = rk;
                    ++vi;
                    if (vi == VT.VecSize)
                    {
                        RealVec vx = VT.same(abc.a.x) * rvi + VT.same(abc.b.x) * rvj + VT.same(abc.c.x) * rvk;
                        RealVec vy = VT.same(abc.a.y) * rvi + VT.same(abc.b.y) * rvj + VT.same(abc.c.y) * rvk;
                        RealVec vz = VT.same(abc.a.z) * rvi + VT.same(abc.b.z) * rvj + VT.same(abc.c.z) * rvk;
                        genericP2L_rowit(make_xyzq(vx, vy, vz, VT.same(1.0)), temp_accu, p);
                        vi = 0;
                    }
                }
            }
        }
    }
    if (vi > 0) {
        RealVec q = VT.same(1.0);
        for ( ; vi < VT.VecSize; ++vi) {
            VT.element(rvi, vi) = VT.element(rvj, vi) = VT.element(rvk, vi) = 1;
            VT.element(q, vi) = 0;
        }
        RealVec vx = VT.same(abc.a.x) * rvi + VT.same(abc.b.x) * rvj + VT.same(abc.c.x) * rvk;
        RealVec vy = VT.same(abc.a.y) * rvi + VT.same(abc.b.y) * rvj + VT.same(abc.c.y) * rvk;
        RealVec vz = VT.same(abc.a.z) * rvi + VT.same(abc.b.z) * rvj + VT.same(abc.c.z) * rvk;
        genericP2L_rowit(make_xyzq(vx, vy, vz, q), temp_accu, p);
    }
    for (ssize_t l = 0; l <= ssize_t(p); ++l)
        for (ssize_t m = 0; m <= l; ++m) {
            Vec_traits<temp_storage_type> VTtemp;
            temp_complex_type cv = accu_vec(l, m);
            accu(l, m) += complex_type(VTtemp.sum(cv.real()), VTtemp.sum(cv.imag()));
        }
    return cnt;
}
#undef ACCVARIANT
#endif  // __AVX__

template <bool enable_x, bool enable_y, bool enable_z,
    typename Real33, typename CoefficientMatrix>
size_t __attribute__((noinline)) L_star_operator(
        const Real33 & abc,
        const size_t ws,
        CoefficientMatrix & l_star,
        const size_t p)
{
    const ssize_t outer = ws * (2 * ws + 1) + ws;
    const ssize_t inner = ws;
//    if (enable_x && enable_y && enable_z)
//        return Lattice3Dcubic_layer_sum<enable_x, enable_y, enable_z>(abc, inner, outer, l_star, p);
#if 0
//#ifdef __AVX__
    return L_star_layer_sum_vec<enable_x, enable_y, enable_z>(abc, inner, outer, l_star, p);
#else
    return L_star_layer_sum_reference<enable_x, enable_y, enable_z>(abc, inner, outer, l_star, p);
#endif
}

}//namespace end

#endif
