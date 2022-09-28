#ifndef _BK_FMM_ostar_operator_hpp_
#define _BK_FMM_ostar_operator_hpp_

#include <functional>
#include "particle2multipole.hpp"
#include "xyzq.hpp"
#include "generic_accumulate.hpp"
#include "super_accumulate.hpp"
#include "iterables.hpp"
#include "triangular_array_adaptor.hpp"
#include "energy_stress_triangular_array_group.hpp"

// for benchmarking and disassembling
#define ostar_cubesum_noinline __attribute__((noinline))

namespace gmx_gpu_fmm{

//
// O_* computation is an *assigning* operation.
// The result will be placed in ostar (a reference to some triangular
// array container), either by copying it there from some temporary
// storage or by first zeroing it and thereafter accumulating into it.
//
template <bool enable_a, bool enable_b, bool enable_c,
    typename Real33, typename CoefficientMatrix>
ostar_cubesum_noinline
size_t O_star_cube_sum_reference(
        const Real33 & abc,
        const size_t ws,
        CoefficientMatrix & ostar,
        const size_t p)
{
    typedef typename Real33::value_type Real3;
    typedef typename Real3::value_type Real;

    assert(p <= ostar.p());
    ostar.zero(p);

    const ssize_t ws_a = enable_a ? ws : 0;
    const ssize_t ws_b = enable_b ? ws : 0;
    const ssize_t ws_c = enable_c ? ws : 0;

    size_t cnt = 0;
    for (ssize_t i = -ws_a; i <= +ws_a; ++i)
        for (ssize_t j = -ws_b; j <= +ws_b; ++j)
            for (ssize_t k = -ws_c; k <= +ws_c; ++k) {
                ++cnt;
                Real3 m = abc * Real3(i, j, k);
                //std::cout<<m<<std::endl;
                P2M(make_xyzq(m, Real(1.0)), ostar, p);  // accumulate
            }
    return cnt;
}


template <bool enable_a, bool enable_b, bool enable_c,
    typename Real33, typename CoefficientMatrix>
ostar_cubesum_noinline
size_t O_star_cube_sum_generic_colstore(
        const Real33 & abc,
        const size_t ws,
        CoefficientMatrix & ostar,
        const size_t p)
{
    typedef typename Real33::value_type Real3;
    typedef typename Real3::value_type Real;
    typedef typename CoefficientMatrix::value_type complex_type;

    typedef simple_triangular_array_adaptor<CoefficientMatrix, complex_type, just_accumulate> storage_adapter_type;

    assert(p <= ostar.p());
    ostar.zero(p);

    const ssize_t ws_a = enable_a ? ws : 0;
    const ssize_t ws_b = enable_b ? ws : 0;
    const ssize_t ws_c = enable_c ? ws : 0;

    P2M_divtable<Real> div_lut_p(p);
    storage_adapter_type storage_adapter(ostar);

    size_t cnt = 0;
    for (ssize_t i = -ws_a; i <= +ws_a; ++i) {
        for (ssize_t j = -ws_b; j <= +ws_b; ++j) {
            for (ssize_t k = -ws_c; k <= +ws_c; ++k) {
                ++cnt;
                Real3 m = abc * Real3(i, j, k);
                P2M_generic(make_xyzq(m, Real(1.0)), storage_adapter, p, div_lut_p);
            }
        }
    }
    return cnt;
}


template <bool enable_a, bool enable_b, bool enable_c,
    typename Real33, typename CoefficientMatrix>
ostar_cubesum_noinline
size_t O_star_cube_sum_generic_rowstore(
        const Real33 & abc,
        const size_t ws,
        CoefficientMatrix & ostar,
        const size_t p)
{
    typedef typename Real33::value_type Real3;
    typedef typename Real3::value_type Real;
    typedef typename CoefficientMatrix::value_type complex_type;

    typedef UpperTriangularArrayRowMajor<complex_type, void> temp_triangle_type;
    typedef simple_triangular_array_adaptor<temp_triangle_type, complex_type, just_accumulate> storage_adapter_type;

    assert(p <= ostar.p());

    const ssize_t ws_a = enable_a ? ws : 0;
    const ssize_t ws_b = enable_b ? ws : 0;
    const ssize_t ws_c = enable_c ? ws : 0;

    P2M_divtable<Real> div_lut_p(p);
    temp_triangle_type temp_storage(p);
    storage_adapter_type storage_adapter(temp_storage);

    size_t cnt = 0;
    for (ssize_t i = -ws_a; i <= +ws_a; ++i) {
        for (ssize_t j = -ws_b; j <= +ws_b; ++j) {
            for (ssize_t k = -ws_c; k <= +ws_c; ++k) {
                ++cnt;
                Real3 m = abc * Real3(i, j, k);
                P2M_generic(make_xyzq(m, Real(1.0)), storage_adapter, p, div_lut_p);
            }
        }
    }
    for (ssize_t l = 0; l <= ssize_t(p); ++l)
        for (ssize_t m = 0; m <= l; ++m)
            ostar(l, m) = temp_storage(l, m);
    return cnt;
}


template <bool enable_a, bool enable_b, bool enable_c,
    typename Real33, typename CoefficientMatrix>
ostar_cubesum_noinline
size_t O_star_cube_sum_rowit(
        const Real33 & abc,
        const size_t ws,
        CoefficientMatrix & ostar,
        const size_t p)
{
    typedef typename Real33::value_type Real3;
    typedef typename Real3::value_type Real;
    typedef typename CoefficientMatrix::value_type complex_type;

    typedef UpperTriangularArrayRowMajor<complex_type, void> temp_triangle_type;
    typedef row_iterable<temp_triangle_type, complex_type, just_accumulate> storage_adapter_type;

    assert(p <= ostar.p());

    const ssize_t ws_a = enable_a ? ws : 0;
    const ssize_t ws_b = enable_b ? ws : 0;
    const ssize_t ws_c = enable_c ? ws : 0;

    P2M_divtable<Real> div_lut_p(p);
    temp_triangle_type temp_storage(p);
    storage_adapter_type storage_adapter(temp_storage);

    size_t cnt = 0;
    for (ssize_t i = -ws_a; i <= +ws_a; ++i) {
        for (ssize_t j = -ws_b; j <= +ws_b; ++j) {
            for (ssize_t k = -ws_c; k <= +ws_c; ++k) {
                ++cnt;
                Real3 m = abc * Real3(i, j, k);
                genericP2M_rowit(make_xyzq(m, Real(1.0)), storage_adapter, p, div_lut_p);
            }
        }
    }
    for (ssize_t l = 0; l <= ssize_t(p); ++l)
        for (ssize_t m = 0; m <= l; ++m)
            ostar(l, m) = temp_storage(l, m);
    return cnt;
}


template <typename Storage, typename Value, template <typename, typename> class Operation>
class scaleaccumulate : public std::binary_function<Storage, Value, void>
{
    typedef Storage storage_type;
    typedef typename storage_type::value_type storage_value_type;
    typedef Value compute_value_type;
    typedef Operation<storage_value_type, compute_value_type> Op;
    typedef typename compute_value_type::value_type scale_type;
public:
    typedef typename storage_type::index_type index_type;

    scaleaccumulate(scale_type scale)
        : scale(scale)
    { }

    void operator () (storage_type & st, index_type l, index_type m, const compute_value_type & value) const
    {
        const Op op;
        op(st(l, m), value * scale);
    }

private:
    scale_type scale;
};


template <bool enable_a, bool enable_b, bool enable_c,
    typename Real33, typename CoefficientMatrix>
ostar_cubesum_noinline
size_t O_star_cube_sum_generic_colstore_scaled(
        const Real33 & abc,
        const size_t ws,
        CoefficientMatrix & ostar,
        const size_t p)
{
    typedef typename Real33::value_type Real3;
    typedef typename Real3::value_type Real;
    typedef typename CoefficientMatrix::value_type complex_type;

    typedef scaleaccumulate<CoefficientMatrix, complex_type, just_accumulate> store_scaled_operator_type;
    typedef custom_triangular_array_adaptor<store_scaled_operator_type> storage_adapter_type;

    assert(p <= ostar.p());

    const ssize_t ws_a = enable_a ? ws : 0;
    const ssize_t ws_b = enable_b ? ws : 0;
    const ssize_t ws_c = enable_c ? ws : 0;

    P2M_divtable<Real> div_lut_p(p);

    size_t cnt = 0;
    for (ssize_t i = -ws_a; i <= +ws_a; ++i) {
        for (ssize_t j = -ws_b; j <= +ws_b; ++j) {
            for (ssize_t k = -ws_c; k <= +ws_c; ++k) {
                ++cnt;
                Real3 m = abc * Real3(i, j, k);
                storage_adapter_type storage_adapter(ostar, Real(k));
                P2M_generic(make_xyzq(m, Real(1.0)), storage_adapter, p, div_lut_p);
            }
        }
    }
    return cnt;
}


template <bool enable_a, bool enable_b, bool enable_c,
    typename Real33, typename ESCoefficientMatrix>
ostar_cubesum_noinline
size_t O_star_es_cube_sum_generic_colstore(
        const Real33 & abc,
        const size_t ws,
        ESCoefficientMatrix & ostar_es,
        const size_t p)
{
    typedef typename Real33::value_type Real3;
    typedef typename Real3::value_type Real;
    typedef typename ESCoefficientMatrix::value_type::value_type complex_type;

    typedef EnergyStressStorageOperator<ESCoefficientMatrix, complex_type, just_accumulate> store_energy_stress_operator_type;
    typedef custom_triangular_array_adaptor<store_energy_stress_operator_type> storage_adapter_type;

    assert(p <= ostar_es.energy.p());
    ostar_es.zero(p);

    const ssize_t ws_a = enable_a ? ws : 0;
    const ssize_t ws_b = enable_b ? ws : 0;
    const ssize_t ws_c = enable_c ? ws : 0;

    P2M_divtable<Real> div_lut_p(p);

    size_t cnt = 0;
    for (ssize_t i = -ws_a; i <= +ws_a; ++i) {
        for (ssize_t j = -ws_b; j <= +ws_b; ++j) {
            for (ssize_t k = -ws_c; k <= +ws_c; ++k) {
                ++cnt;
                Real3 ijk(i, j, k);
                Real3 m = abc * -ijk;
                storage_adapter_type storage_adapter(ostar_es, ijk.x, ijk.y, ijk.z);
                P2M_generic(make_xyzq(m, Real(1.0)), storage_adapter, p, div_lut_p);
            }
        }
    }
    return cnt;
}


template <bool enable_a, bool enable_b, bool enable_c,
    typename Real33, typename ESCoefficientMatrix>
ostar_cubesum_noinline
size_t O_star_es_cube_sum_generic_rowstore(
        const Real33 & abc,
        const size_t ws,
        ESCoefficientMatrix & ostar_es,
        const size_t p)
{
    typedef typename Real33::value_type Real3;
    typedef typename Real3::value_type Real;
    typedef typename ESCoefficientMatrix::value_type::value_type complex_type;

    typedef UpperTriangularArrayRowMajor<complex_type, void> temp_triangular_type;
    typedef EnergyStressTriangularArrayGroup<temp_triangular_type> es_temp_triangular_type;
    typedef EnergyStressStorageOperator<es_temp_triangular_type, complex_type, just_accumulate> store_energy_stress_operator_type;
    typedef custom_triangular_array_adaptor<store_energy_stress_operator_type> storage_adapter_type;

    assert(p <= ostar_es.energy.p());

    const ssize_t ws_a = enable_a ? ws : 0;
    const ssize_t ws_b = enable_b ? ws : 0;
    const ssize_t ws_c = enable_c ? ws : 0;

    P2M_divtable<Real> div_lut_p(p);
    es_temp_triangular_type temp_storage(p);

    size_t cnt = 0;
    for (ssize_t i = -ws_a; i <= +ws_a; ++i) {
        for (ssize_t j = -ws_b; j <= +ws_b; ++j) {
            for (ssize_t k = -ws_c; k <= +ws_c; ++k) {
                ++cnt;
                Real3 ijk(i, j, k);
                Real3 m = abc * -ijk;
                storage_adapter_type storage_adapter(temp_storage, ijk.x, ijk.y, ijk.z);
                P2M_generic(make_xyzq(m, Real(1.0)), storage_adapter, p, div_lut_p);
            }
        }
    }
    for (ssize_t l = 0; l <= ssize_t(p); ++l)
        for (ssize_t m = 0; m <= l; ++m) {
            ostar_es.energy(l, m) = temp_storage.energy(l, m);
            ostar_es.stress.a(l, m) = temp_storage.stress.a(l, m);
            ostar_es.stress.b(l, m) = temp_storage.stress.b(l, m);
            ostar_es.stress.c(l, m) = temp_storage.stress.c(l, m);
        }
    return cnt;
}

#if 0
//#ifdef __AVX__
template <bool enable_x, bool enable_y, bool enable_z,
    typename Real33, typename CoefficientMatrix>
ostar_cubesum_noinline
size_t O_star_cube_sum_vec(
        const Real33 & abc,
        const size_t ws,
        CoefficientMatrix & ostar,
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
#elif ACCVARIANT == 2
    // compute with SuperSIMD, temporary storage: SIMD duplicated interleaved IRIR
    typedef simd_type temp_storage_type;
#else
    // compute with SuperSIMD, temporary storage: scalar
    typedef Real temp_storage_type;
#endif
    typedef typename COMPLEX_GENERATOR<RealVec>::type simd_complex_type;
#if ACCVARIANT == 2
    typedef temp_storage_type temp_complex_type;
#else
    typedef typename COMPLEX_GENERATOR<temp_storage_type>::type temp_complex_type;
#endif
    typedef aligned_allocator<temp_complex_type> temp_allocator_type;
    typedef UpperTriangularArrayRowMajor<temp_complex_type, void, temp_allocator_type> temp_triangle_type;
#if ACCVARIANT == 0
    typedef row_iterable<temp_triangle_type, simd_complex_type, just_accumulate> temp_accumulator_type;
#elif ACCVARIANT == 1
    typedef row_iterable<temp_triangle_type, simd_complex_type, just_complex_super_accumulate> temp_accumulator_type;
#elif ACCVARIANT == 2
    typedef row_iterable<temp_triangle_type, simd_complex_type, just_complex_interleaved_super_accumulate> temp_accumulator_type;
#else
    typedef row_iterable<temp_triangle_type, simd_complex_type, just_complex_scalar_accumulate> temp_accumulator_type;
#endif

    assert(p <= ostar.p());

    const ssize_t ws_x = enable_x ? ws : 0;
    const ssize_t ws_y = enable_y ? ws : 0;
    const ssize_t ws_z = enable_z ? ws : 0;

    Vec_traits<RealVec> VT;
    RealVec rvi, rvj, rvk;
    size_t vi = 0;
    temp_triangle_type accu_vec(p);
    temp_accumulator_type temp_accu(accu_vec);
#if 1
    P2M_divtable<RealVec> div_lut_p(p);
#else
    P2M_divtable<RealVec, true, aligned_allocator<RealVec> > div_lut_p(p);
#endif

    size_t cnt = 0;
    for (ssize_t i = -ws_x; i <= +ws_x; ++i) {
        Real ri(i);
        for (ssize_t j = -ws_y; j <= +ws_y; ++j) {
            Real rj(j);
            for (ssize_t k = -ws_z; k <= +ws_z; ++k) {
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
                    genericP2M_rowit(make_xyzq(vx, vy, vz, VT.same(1.0)), temp_accu, p, div_lut_p);
                    vi = 0;
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
        genericP2M_rowit(make_xyzq(vx, vy, vz, q), temp_accu, p, div_lut_p);
    }
    for (ssize_t l = 0; l <= ssize_t(p); ++l)
        for (ssize_t m = 0; m <= l; ++m) {
            Vec_traits<temp_storage_type> VTtemp;
            temp_complex_type cv = accu_vec(l, m);
#if ACCVARIANT == 2
// TODO
#else
            ostar(l, m) = complex_type(VTtemp.sum(cv.real()), VTtemp.sum(cv.imag()));
#endif
        }
    return cnt;
}
#undef ACCVARIANT
#endif  // __AVX__

template <bool enable_x, bool enable_y, bool enable_z,
    typename Real33, typename CoefficientMatrix>
__attribute__((noinline))
size_t O_star_operator(
        const Real33 & abc,
        const size_t ws,
        CoefficientMatrix & o_star,
        const size_t p)
{
#if 0
//#ifdef __AVX__
    return O_star_cube_sum_vec<enable_x, enable_y, enable_z>(abc, ws, o_star, p);
#else
    return O_star_cube_sum_reference<enable_x, enable_y, enable_z>(abc, ws, o_star, p);
#endif
}

}//namespace end

#endif
// vim: et:ts=4:sw=4
