#ifndef _BK_FMM_multipole_hpp_
#define _BK_FMM_multipole_hpp_

#include "fmm_complex.hpp"
#include "triangular_array.hpp"
#include "parity_sign.hpp"
#include "cuda_alloc.hpp"
#include "architecture.hpp"

namespace gmx_gpu_fmm{

template <typename Real, typename arch>
class MultipoleCoefficientsUpper : public UpperTriangularArray<typename COMPLEX_GENERATOR<Real>::type, arch>
{
    typedef UpperTriangularArray<typename COMPLEX_GENERATOR<Real>::type, arch> Base;
public:
    typedef arch architecture;
    typedef Real real_type;
    typedef typename Base::value_type complex_type;
    typedef typename Base::index_type index_type;


    MultipoleCoefficientsUpper(index_type p) : Base(p)
    {}

    MultipoleCoefficientsUpper()
    {}
};

template <typename Real, typename arch>
class MultipoleCoefficientsUpperSoA : public UpperTriangularArray<typename COMPLEX_GENERATOR<Real>::pointer, arch>
{
    typedef typename COMPLEX_GENERATOR<Real>::pointer pointer_to_data_type;
    typedef typename COMPLEX_GENERATOR<Real>::type data_type;

    typedef UpperTriangularArray<pointer_to_data_type, arch> Base;
public:
    typedef Real real_type;
    typedef typename Base::value_type complex_type;
    typedef typename Base::index_type index_type;


    MultipoleCoefficientsUpperSoA(index_type p) : Base(p)
    {}

    MultipoleCoefficientsUpperSoA()
    {}
};

template <typename Real, typename arch>
class MultipoleCoefficientsUpper_Memory : public UpperTriangularArray_Memory<typename COMPLEX_GENERATOR<Real>::type, arch>
{
    typedef UpperTriangularArray_Memory<typename COMPLEX_GENERATOR<Real>::type, arch> Base;
public:
    typedef Real real_type;
    typedef typename Base::value_type complex_type;
    typedef typename Base::index_type index_type;


    MultipoleCoefficientsUpper_Memory(index_type p, size_t num_of) : Base(p, num_of)
    { }

    MultipoleCoefficientsUpper_Memory(){}
};

template <typename Real>
class MultipoleCoefficientsUpperLower : public TriangularArray<typename COMPLEX_GENERATOR<Real>::type, cuda_allocator<typename COMPLEX_GENERATOR<Real>::type> >
{
    typedef TriangularArray<typename COMPLEX_GENERATOR<Real>::type, cuda_allocator<typename COMPLEX_GENERATOR<Real>::type> > Base;
public:
    typedef Real real_type;
    typedef typename Base::value_type complex_type;
    typedef typename Base::index_type index_type;

    MultipoleCoefficientsUpperLower(index_type p) : Base(p)
    {}

    template <typename Table>

    MultipoleCoefficientsUpperLower(index_type p, const Table & table) : Base(p)
    {
        for (size_t i = 0; i < sizeof(table) / sizeof(table[0]); ++i)
        {
            index_type l = table[i].l, m = table[i].m;
            if (l >= 0 && l <= p && m >= -l && m <= l)
                this->get(l, m) = table[i].val;
        }
        this->populate_lower();
    }
};

}//namespace end

#endif
