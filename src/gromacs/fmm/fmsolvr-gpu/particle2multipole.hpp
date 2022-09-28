#ifndef _BK_FMM_particle2multipole_hpp_
#define _BK_FMM_particle2multipole_hpp_

#include "Vec.h"
#include "global_functions.hpp"
#include "xyzq.hpp"

namespace gmx_gpu_fmm{

template <typename StorageT, typename ValueT, bool store_simd>
struct divtable_load_expand;

template <typename StorageT, typename ValueT>
struct divtable_load_expand<StorageT, ValueT, true>
{
    ValueT operator () (const StorageT & val) const
    {
        return val;
    }
};

template <typename StorageT, typename ValueT>
struct divtable_load_expand<StorageT, ValueT, false>
{
    ValueT operator () (const StorageT & val) const
    {
        Vec_traits<ValueT> VT;
        return VT.same(val);
    }
};

template <typename StorageT, typename ValueT, typename Expander>
class divtable_const_iterator
{
    typedef divtable_const_iterator self;
    typedef StorageT storage_type;
public:
    typedef ValueT value_type;

    divtable_const_iterator(const storage_type * p) : p(p)
    { }

    self & operator ++ ()
    {
        ++p;
        return *this;
    }

    self operator ++ (int)
    {
        self ret(*this);
        ++p;
        return ret;
    }

    value_type operator * () const
    {
        return Expander()(*p);
    }

private:
    const storage_type * p;
};

template<typename T, typename U, bool store_simd>
struct conditional
{
  typedef typename T::type storage_type;
};

template<typename T, typename U>
struct conditional<T,U,true>
{
  typedef typename U::type storage_type;
};

template <typename Real, bool store_simd = false, typename Allocator = std::allocator<Real> >
class P2M_divtable {
    typedef Vec_traits<Real> realvectraits;
    typedef typename realvectraits::scalar_type scalar_type;
    typedef typename realvectraits::vec_type vec_type;
    
    //typedef typename std::conditional<store_simd, vec_type, scalar_type>::type storage_type;

    typedef typename conditional<scalar_type,vec_type,store_simd>::storage_type storage_type;
      
    typedef typename Allocator::template rebind<storage_type>::other allocator_type;
    typedef size_t size_type;
    typedef divtable_load_expand<storage_type, vec_type, store_simd> load_expander;
public:
    typedef divtable_const_iterator<storage_type, vec_type, load_expander> const_iterator;

    //__attribute__((noinline))
    P2M_divtable(size_type p)
        : p_(p)
    {
        table = allocator.allocate(size());
        for (size_type i = 0; i < size(); ++i)
            allocator.construct(&table[i], storage_type());

        Vec_traits<scalar_type> ST;
        size_type i = 0;

        for (size_type m = 0; m < p; ++m) {
            store(i++, reciprocal(ST.same(2 * (m + 1))));

            for (size_type l = m + 2; l <= p; ++l) {
                store(i++, reciprocal(ST.same(l * l - m * m)));
            }
        }
    }

    ~P2M_divtable()
    {
        for (size_type i = 0; i < size(); ++i)
            allocator.destroy(&table[i]);
        allocator.deallocate(table, size());
    }

    size_t p() const
    {
        return p_;
    }

    const_iterator begin() const
    {
        return const_iterator(table);
    }

#if 1
    vec_type operator [] (size_type i) const
    {
        return load_expander()(table[i]);
        //return expand(table[i]);
    }
#else
    template <bool StoreSimd = store_simd>
    typename std::enable_if<!StoreSimd, vec_type>::type
        operator [] (size_type i) const
    {
        realvectraits VT;
        return VT.same(table[i]);
    }

    template <bool StoreSimd = store_simd>
    typename std::enable_if<StoreSimd, vec_type>::type
        operator [] (size_type i) const
    {
        return table[i];
    }
#endif

private:
    storage_type * table;
    allocator_type allocator;
    size_t p_;

    size_t size() const
    {
        return p_ * (p_ + 1) / 2;
    }

    /*
    template <bool StoreSimd = store_simd>
    typename std::enable_if<!StoreSimd, void>::type
        store(size_type i, scalar_type val)
    {
        table[i] = val;
    }
    */
    
    //template <bool StoreSimd = store_simd>
    //typename std::enable_if<StoreSimd, void>::type
    void
        store(size_type i, scalar_type val)
    {
        realvectraits VT;
        table[i] = VT.same(val);
    }
    /*
    template <bool StoreSimd = store_simd>
    typename std::enable_if<!StoreSimd, vec_type>::type
        expand(const storage_type & val) const
    {
        realvectraits VT;
        return VT.same(val);
    }
    */
    
    //template <bool StoreSimd = store_simd>
    //typename std::enable_if<StoreSimd, vec_type>::type
    vec_type 
       expand(const storage_type & val) const
        
    {
        return val;
    }
   
};

// add multipole coefficients for a particle with charge q at relative position (x,y,z)
// to the (non-negative triangle) multipole coefficient matrix upto order p
template <typename Real4, typename CoefficientMatrix>
inline void P2M_old(
        const Real4 & xyzq,
        CoefficientMatrix & omega,
        const size_t p)
{
    printf("running old p2m?\n");
    typedef typename Real4::value_type real_type;
    typedef typename CoefficientMatrix::value_type complex_type;
    Vec_traits<real_type> VT;

    const real_type x = xyzq.x;
    const real_type y = xyzq.y;
    const real_type z = xyzq.z;
    const real_type q = xyzq.q;

    real_type dist_squared = x * x + y * y + z * z;

    // omega_0_0  (for the first iteration)
    complex_type omega_mplus1_mplus1(q);

    // omega_0_0 upto omega_p_p-1
    for (size_t m = 0; m < p; ++m) {
        // omega_m_m  (from previous iteration)
        complex_type omega_m_m = omega_mplus1_mplus1;
        omega(m, m) += omega_m_m;

        // omega_m+1_m+1  (for the next iteration)
        omega_mplus1_mplus1 = reciprocal(VT.same(2 * (m + 1))) * complex_type(x, -y) * omega_m_m;

        // omega_m+1_m
        complex_type omega_mplus1_m = z * omega_m_m;
        omega(m + 1, m) += omega_mplus1_m;

        // omega_m+2_m upto omega_p_m
        complex_type omega_lminus2_m = omega_m_m;
        complex_type omega_lminus1_m = omega_mplus1_m;
        for (size_t l = m + 2; l <= p; ++l) {
            // omega_l_m
            complex_type omega_l_m = reciprocal(VT.same(l * l - m * m)) *
                (VT.same(2 * l - 1) * z * omega_lminus1_m - dist_squared * omega_lminus2_m);
            omega(l, m) += omega_l_m;
            omega_lminus2_m = omega_lminus1_m;
            omega_lminus1_m = omega_l_m;
        }
    }

    // omega_p_p
    omega(p, p) += omega_mplus1_mplus1;
}


// add multipole coefficients for a particle with charge q at relative position (x,y,z)
// to the (non-negative triangle) multipole coefficient matrix upto order p
template <typename Real4, typename CoefficientMatrix>
inline void P2M_delayed_complex(
        const Real4 & xyzq,
        CoefficientMatrix & omega,
        const size_t p)
{

    printf("running delayed p2m?\n");
    typedef typename Real4::value_type real_type;
    typedef typename CoefficientMatrix::value_type complex_type;
    typedef typename Vec_traits<real_type>::scalar_type scalar_type;
    Vec_traits<real_type> VT;

    const real_type x = xyzq.x;
    const real_type y = xyzq.y;
    const real_type z = xyzq.z;
    const real_type q = xyzq.q;

    real_type dist_squared = x * x + y * y + z * z;

    complex_type complex_x_y(x, -y);
    complex_type complex_x_y_m(VT.same(1));     // iterative computation of ((x - iy)^m)

    // omega_0_0  (for the first iteration)
    real_type omega_mplus1_mplus1 = q;

    // omega_0_0 upto omega_p_p-1
    for (size_t m = 0; m < p; ++m) {
        // omega_m_m  (from previous iteration)
        real_type omega_m_m = omega_mplus1_mplus1;
        omega(m, m) += complex_x_y_m * omega_m_m;

        // omega_m+1_m+1  (for the next iteration)
        omega_mplus1_mplus1 = VT.same(reciprocal(scalar_type(2 * (m + 1)))) * omega_m_m;

        // omega_m+1_m
        real_type omega_mplus1_m = z * omega_m_m;
        omega(m + 1, m) += complex_x_y_m * omega_mplus1_m;

        // omega_m+2_m upto omega_p_m
        real_type omega_lminus2_m = omega_m_m;
        real_type omega_lminus1_m = omega_mplus1_m;
        for (size_t l = m + 2; l <= p; ++l) {
            // omega_l_m
            real_type omega_l_m = VT.same(reciprocal(scalar_type(l * l - m * m))) *
                (VT.same(2 * l - 1) * z * omega_lminus1_m - dist_squared * omega_lminus2_m);
            omega(l, m) += complex_x_y_m * omega_l_m;

            // (for the next iteration)
            omega_lminus2_m = omega_lminus1_m;
            omega_lminus1_m = omega_l_m;
        }

        // (for the next iteration)
        complex_x_y_m *= complex_x_y;
    }

    // omega_p_p
    omega(p, p) += complex_x_y_m * omega_mplus1_mplus1;
}


// add multipole coefficients for a particle with charge q at relative position (x,y,z)
// to the (non-negative triangle) multipole coefficient matrix upto order p
template <typename Real4, typename CoefficientMatrix>
inline void P2M(
        const Real4 & xyzq,
        CoefficientMatrix & omega,
        const size_t p)
{
    P2M_nodivtable(xyzq.q, xyzq, omega, p);
}

// add multipole coefficients for a particle with charge q at relative position (x,y,z)
// to the (non-negative triangle) multipole coefficient matrix upto order p

template <typename T>
T set_min_max2(T val)
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

template <typename Real3, typename CoefficientMatrix>
__attribute__((noinline))
__attribute__((flatten))
void P2M_nodivtable(
        const typename Real3::value_type & q_,
        const Real3 & xyz,
        CoefficientMatrix & omega,
        const size_t p)
{
    //printf("running nodiv p2m?\n");
    //typedef double real_type;
    //typedef fmm::complex<real_type> complex_type;
    typedef typename Real3::value_type real_type;
    typedef typename CoefficientMatrix::value_type complex_type;
    typedef typename Vec_traits<real_type>::scalar_type scalar_type;
    Vec_traits<real_type> VT;

    const real_type x = xyz.x;
    const real_type y = xyz.y;
    const real_type z = xyz.z;
    const real_type q = q_;

    real_type dist_squared = x * x + y * y + z * z;
    real_type twice_z = z + z;

    complex_type complex_x_y(x, -y);
    complex_type complex_x_y_m(VT.same(1));     // iterative computation of ((x - iy)^m)
    

    // omega_0_0  (for the first iteration)
    real_type omega_mplus1_mplus1 = q;
    real_type e_m = z + twice_z;                // iterative computation of ((2*m + 3)*z)

    // omega_0_0 upto omega_p_p-1
    for (size_t m = 0; m < p; ++m) {
        // omega_m_m  (from previous iteration)
        real_type omega_m_m = omega_mplus1_mplus1;
        omega(m, m) += set_min_max2(complex_x_y_m * omega_m_m);

        // omega_m+1_m+1  (for the next iteration)
        omega_mplus1_mplus1 = VT.same(reciprocal(scalar_type(2 * (m + 1)))) * omega_m_m;

        // omega_m+1_m
        real_type omega_mplus1_m = z * omega_m_m;
        omega(m + 1, m) += set_min_max2(complex_x_y_m * omega_mplus1_m);

        // omega_m+2_m upto omega_p_m
        real_type omega_lminus2_m = omega_m_m;
        real_type omega_lminus1_m = omega_mplus1_m;
        real_type f_l = e_m;                    // iterative computation of ((2*l - 1)*z)
        for (size_t l = m + 2; l <= p; ++l) {
            // omega_l_m
            real_type omega_l_m = VT.same(reciprocal(scalar_type(l * l - m * m))) *
                (f_l * omega_lminus1_m - dist_squared * omega_lminus2_m);
            omega(l, m) += set_min_max2(complex_x_y_m * omega_l_m);

            // (for the next iteration)
            omega_lminus2_m = omega_lminus1_m;
            omega_lminus1_m = omega_l_m;
            f_l += twice_z;
        }

        // (for the next iteration)
        complex_x_y_m *= complex_x_y;
        e_m += twice_z;
    }

    // omega_p_p
    omega(p, p) += set_min_max2(complex_x_y_m * omega_mplus1_mplus1);
}


// add multipole coefficients for a particle with charge q at relative position (x,y,z)
// to the (non-negative triangle) multipole coefficient matrix upto order p
template <typename Real3, typename CoefficientMatrix, typename DivisionTable>
inline
__attribute__((flatten))
//__attribute__((noinline))
void P2M_reference(
        const typename Real3::value_type & q,
        const Real3 & xyz,
        CoefficientMatrix & omega,
        const DivisionTable & div_lut)
{
    printf("running reference p2m?\n");
    typedef typename Real3::value_type real_type;
    typedef typename CoefficientMatrix::value_type complex_type;

    Vec_traits<real_type> VT;
    const size_t p = div_lut.p();
    assert(p <= omega.p());

    const real_type x = xyz.x;
    const real_type y = xyz.y;
    const real_type z = xyz.z;

    real_type dist_squared = x * x + y * y + z * z;
    real_type twice_z = z + z;
    size_t i = 0;

    complex_type complex_x_y(x, -y);
    complex_type complex_x_y_m(VT.same(1));     // iterative computation of ((x - iy)^m)

    // omega_0_0  (for the first iteration)
    real_type omega_mplus1_mplus1 = q;
    real_type e_m = z + twice_z;                // iterative computation of ((2*m + 3)*z)

    // omega_0_0 upto omega_p_p-1
    for (size_t m = 0; m < p; ++m) {
        // omega_m_m  (from previous iteration)
        real_type omega_m_m = omega_mplus1_mplus1;
        omega(m, m) += complex_x_y_m * omega_m_m;

        // omega_m+1_m+1  (for the next iteration)
        omega_mplus1_mplus1 = div_lut[i++] * omega_m_m;

        // omega_m+1_m
        real_type omega_mplus1_m = z * omega_m_m;
        omega(m + 1, m) += complex_x_y_m * omega_mplus1_m;

        // omega_m+2_m upto omega_p_m
        real_type omega_lminus2_m = omega_m_m;
        real_type omega_lminus1_m = omega_mplus1_m;
        real_type f_l = e_m;                    // iterative computation of ((2*l - 1)*z)


        for (size_t l = m + 2; l <= p; ++l) {
            printf("%lu ",l);
            // omega_l_m
            real_type omega_l_m = div_lut[i++] *
                (f_l * omega_lminus1_m - dist_squared * omega_lminus2_m);
            omega(l, m) += complex_x_y_m * omega_l_m;

            // (for the next iteration)
            omega_lminus2_m = omega_lminus1_m;
            omega_lminus1_m = omega_l_m;
            f_l += twice_z;
        }
        printf("M %lu \n",m);

        // (for the next iteration)
        complex_x_y_m *= complex_x_y;
        e_m += twice_z;
    }

    // omega_p_p
    omega(p, p) += complex_x_y_m * omega_mplus1_mplus1;
}


// add multipole coefficients for a particle with charge q at relative position (x,y,z)
// to the (non-negative triangle) multipole coefficient matrix upto order p
template <typename Real4, typename CoefficientMatrix, typename DivisionTable>
inline
void P2M(
        const Real4 & xyzq,
        CoefficientMatrix & omega,
        const size_t p,
        const DivisionTable & div_lut_p)
{
    P2M_reference(xyzq.q, xyzq, omega, div_lut_p);
}


// add multipole coefficients upto order p for the particle at
// 'particle_coords' with strength (charge/mass/...) 'strength' to the
// (non-negative triangle) multipole coefficient matrix
// for a multipole expansion 'omega' at 'expansion_point'
template <typename Real3, typename CoefficientMatrix, typename DivisionTable>
inline
void P2M(
        const typename Real3::value_type & strength,
        const Real3 & particle_coords,
        CoefficientMatrix & omega,
        const Real3 & expansion_point,
        const DivisionTable & div_lut_p)
{
    P2M_reference(strength, particle_coords - expansion_point, omega, div_lut_p);
}


// compute multipole coefficients upto order p for a particle with charge q
// at relative position (x,y,z) and store (accumulate/assign/...) them in
// a (non-negative triangle) multipole coefficient matrix
template <typename Real4, typename CoefficientStorage, typename DivisionTable>
inline
__attribute__((flatten))
//__attribute__((noinline))
void genericP2M_reference(
        const Real4 & xyzq,
        CoefficientStorage & omega,
        const size_t p,
        const DivisionTable & div_lut)
{
    printf("running generic p2m?\n");
    typedef typename CoefficientStorage::value_type complex_type;
    typedef typename complex_type::value_type real_type;
    typedef typename DivisionTable::const_iterator div_table_iterator;
    Vec_traits<real_type> VT;
    assert(p == div_lut.p());
    div_table_iterator div_table = div_lut.begin();

    const real_type x = xyzq.x;
    const real_type y = xyzq.y;
    const real_type z = xyzq.z;
    const real_type q = xyzq.q;

    real_type dist_squared = x * x + y * y + z * z;
    real_type twice_z = z + z;

    complex_type complex_x_y(x, -y);
    complex_type complex_x_y_m(VT.same(1));     // iterative computation of ((x - iy)^m)

    // omega_0_0  (for the first iteration)
    real_type omega_mplus1_mplus1 = q;
    real_type e_m = z + twice_z;                // iterative computation of ((2*m + 3)*z)

    // omega_0_0 upto omega_p_p-1
    for (size_t m = 0; m < p; ++m) {
        // omega_m_m  (from previous iteration)
        real_type omega_m_m = omega_mplus1_mplus1;
        omega(m, m, complex_x_y_m * omega_m_m);

        // omega_m+1_m+1  (for the next iteration)
        omega_mplus1_mplus1 = (*(div_table++) /*  1 / (2 * (m + 1))  */) * omega_m_m;

        // omega_m+1_m
        real_type omega_mplus1_m = z * omega_m_m;
        omega(m + 1, m, complex_x_y_m * omega_mplus1_m);

        // omega_m+2_m upto omega_p_m
        real_type omega_lminus2_m = omega_m_m;
        real_type omega_lminus1_m = omega_mplus1_m;
        real_type f_l = e_m;                    // iterative computation of ((2*l - 1)*z)
        for (size_t l = m + 2; l <= p; ++l) {
            // omega_l_m
            real_type omega_l_m = (*(div_table++) /*  1 / (l * l - m * m)  */) *
                (f_l * omega_lminus1_m - dist_squared * omega_lminus2_m);
            omega(l, m, complex_x_y_m * omega_l_m);

            // (for the next iteration)
            omega_lminus2_m = omega_lminus1_m;
            omega_lminus1_m = omega_l_m;
            f_l += twice_z;
        }

        // (for the next iteration)
        complex_x_y_m *= complex_x_y;
        e_m += twice_z;
    }

    // omega_p_p
    omega(p, p, complex_x_y_m * omega_mplus1_mplus1);
}


// compute multipole coefficients upto order p for a particle with charge q
// at relative position (x,y,z) and store (accumulate/assign/...) them in
// a (non-negative triangle) multipole coefficient matrix
template <typename Real4, typename CoefficientStorage, typename DivisionTable>
inline
__attribute__((flatten))
//__attribute__((noinline))
void genericP2M_rowit(
        const Real4 & xyzq,
        CoefficientStorage & omega,
        const size_t p,
        const DivisionTable & div_lut)
{
    printf("running generic rowit p2m?\n");
    typedef typename CoefficientStorage::row_iterator row_iterator;
    typedef typename CoefficientStorage::value_type complex_type;
    typedef typename complex_type::value_type real_type;
    typedef typename DivisionTable::const_iterator div_table_iterator;
    Vec_traits<real_type> VT;
    assert(p == div_lut.p());
    div_table_iterator div_table = div_lut.begin();

    const real_type x = xyzq.x;
    const real_type y = xyzq.y;
    const real_type z = xyzq.z;
    const real_type q = xyzq.q;

    real_type dist_squared = x * x + y * y + z * z;
    real_type twice_z = z + z;

    complex_type complex_x_y(x, -y);
    complex_type complex_x_y_m(VT.same(1));     // iterative computation of ((x - iy)^m)

    // omega_0_0  (for the first iteration)
    real_type omega_mplus1_mplus1 = q;
    real_type e_m = z + twice_z;                // iterative computation of ((2*m + 3)*z)

    // omega_0_0 upto omega_p_p-1
    for (size_t m = 0; m < p; ++m) {
        row_iterator omega_row_m = omega.get_row_iterator(m, m);

        // omega_m_m  (from previous iteration)
        real_type omega_m_m = omega_mplus1_mplus1;
        (*omega_row_m)(complex_x_y_m * omega_m_m);
        ++omega_row_m;

        // omega_m+1_m+1  (for the next iteration)
        omega_mplus1_mplus1 = (*div_table) * omega_m_m;
        ++div_table;

        // omega_m+1_m
        real_type omega_mplus1_m = z * omega_m_m;
        (*omega_row_m)(complex_x_y_m * omega_mplus1_m);
        ++omega_row_m;

        // omega_m+2_m upto omega_p_m
        real_type omega_lminus2_m = omega_m_m;
        real_type omega_lminus1_m = omega_mplus1_m;
        real_type f_l = e_m;                    // iterative computation of ((2*l - 1)*z)
        for (size_t l = m + 2; l <= p; ++l) {
            // omega_l_m
            real_type omega_l_m = (*div_table) *
                (f_l * omega_lminus1_m - dist_squared * omega_lminus2_m);
            ++div_table;
            (*omega_row_m)(complex_x_y_m * omega_l_m);
            ++omega_row_m;

            // (for the next iteration)
            omega_lminus2_m = omega_lminus1_m;
            omega_lminus1_m = omega_l_m;
            f_l += twice_z;
        }

        // (for the next iteration)
        complex_x_y_m *= complex_x_y;
        e_m += twice_z;
    }

    // omega_p_p
    row_iterator omega_p_p = omega.get_row_iterator(p, p);
    (*omega_p_p)(complex_x_y_m * omega_mplus1_mplus1);
}


// DEPRECATED name
template <typename Real4, typename CoefficientStorage, typename DivisionTable>
inline
void P2M_generic(
        const Real4 & xyzq,
        CoefficientStorage & omega,
        const size_t p,
        const DivisionTable & div_lut)
{
    genericP2M_reference(xyzq, omega, p, div_lut);
}


// Particle2Multipole
// add coefficients contributed by the given particles to omega
// (input omega must be either zeroed or contain proper coefficients to be added to)
template <typename TableReal, typename Real3, typename CoefficientMatrix>
void P2M(
        const TableReal * vx,
        const TableReal * vy,
        const TableReal * vz,
        const TableReal * vq,
        const Real3 & expansion_point,
        const size_t begin,
        const size_t end,
        CoefficientMatrix & omega,
        const size_t p)
{
    typedef typename Real3::value_type Real;

    P2M_divtable<Real> div_lut_p(p);
    for (size_t particle = begin; particle < end; ++particle) {
        Real3 xyz(vx[particle], vy[particle], vz[particle]);
        Real q = vq[particle];
        P2M(q, xyz, omega, expansion_point, div_lut_p);
    }
}

}//namespace end

#endif

