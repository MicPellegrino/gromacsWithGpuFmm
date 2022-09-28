#ifndef _BK_FMM_triangular_array_hpp_
#define _BK_FMM_triangular_array_hpp_

#include <cassert>
#include "cuda_keywords.hpp"
#include <ios>
#include <iomanip>
#include "data_type.hpp"
#include "managed.hpp"
#include "parity_sign.hpp"
#include "architecture.hpp"
#include <limits>

namespace gmx_gpu_fmm{

template <typename T>
struct remove_ptr
{
    typedef T type;
};

template<typename T>
struct remove_ptr<T*>
{
    typedef typename remove_ptr<T>::type type;
};

template<typename T>
struct is_ptr
{
    static const bool value = false;
};

template<typename T>
struct is_ptr<T*>
{
    static const bool value = true;
};

template <typename T, typename architecture>
class UpperTriangularArray: public Managed<architecture>
{

public:

    typedef T value_type;
    typedef typename remove_ptr<T>::type complex_type;
    typedef size_t size_type;
    typedef ssize_t index_type;
    typedef typename complex_type::value_type Realtype;
    typedef typename architecture::template rebind<value_type>::other allocator_type;

    UpperTriangularArray(){}

    UpperTriangularArray(index_type p):p_(p)
    {
        array = allocator.allocate(size());
        for (size_type i = 0; i < size(); ++i)
            allocator.construct(&array[i], value_type());
    }

    void init(index_type p){
        p_ = p;
        array = allocator.allocate(size());
        for (size_type i = 0; i < size(); ++i)
            allocator.construct(&array[i], value_type());
    }

    void reinit(index_type p, value_type *ptr){
        p_= p;
        array = ptr;
    }

    ~UpperTriangularArray()
    {
        for (size_type i = 0; i < size(); ++i)
            allocator.destroy(&array[i]);
        allocator.deallocate(array, size());
    }

    CUDA
    size_type p() const
    {
        return p_;
    }

    CUDA
    void fill(const value_type & value = value_type())
    {
        for(index_type i=0;i<size(); ++i)
        {
            array[i] = value;
        }
    }

    CUDA
    void flush(index_type i)
    {
        array[i] = value_type();
    }

    CUDA
    void zero()
    {
        fill();
    }

    CUDA
    void zero(index_type p)
    {

        for(index_type i=0; i < size(p); ++i)
        {
            array[i] = value_type();
        }
    }

    CUDA
    void scale(Realtype scaling_factor)
    {

        for(index_type i=0; i < size(p_); ++i)
        {
            array[i] *= scaling_factor;
        }
    }

    Realtype max()
    {
        Realtype m = std::abs(std::numeric_limits<Realtype>::min());
        for(index_type i=0; i < size(); ++i)
        {
            if(std::abs(array[i].real()) > m)
                m = std::abs(array[i].real());

            if(std::abs(array[i].imag()) > m)
                m = std::abs(array[i].imag());
        }
        return m;
    }

    Realtype min()
    {
        Realtype m = std::abs(std::numeric_limits<Realtype>::max());
        for(index_type i=0; i < size(); ++i)
        {
            if(std::abs(array[i].real() == 0.0))
                continue;
            if(std::abs(array[i].real()) < m)
                m = std::abs(array[i].real());
            if(std::abs(array[i].imag() == 0.0))
                continue;
            if(std::abs(array[i].imag()) < m)
                m = std::abs(array[i].imag());
        }
        return m;
    }

    // 0 <= l <= p, 0 <= m <= l
    CUDA
    const value_type& get_upper(index_type l, index_type m) const
    {
        return array[offset(l, m)];
    }

    // 0 <= l <= p, 0 <= m <= l
    CUDA
    value_type & get_upper(index_type l, index_type m)
    {
        return array[offset(l, m)];
    }

    // 0 <= l <= p, -l <= m <= -1
    CUDA
    value_type get_lower(index_type l, index_type m) const
    {
        return toggle_sign_if_odd(-m, conj(get_upper(l, -m)));
    }

    // 0 <= l <= p, -l <= m <= l
    CUDA
    value_type get(index_type l, index_type m) const
    {
        return (m < 0) ? get_lower(l, m) : get_upper(l, m);
    }

    CUDA
    value_type& getSoA(index_type l, index_type m)
    {
        return array[offset(l, m)];
    }

    CUDA
    complex_type getSoA(index_type l, index_type m, size_t box_id)
    {

        const complex_type* ptr = array[offset(l, m)];
        const REAL2* tmp_p = reinterpret_cast<const REAL2*>(ptr + box_id);
        const REAL2 tmp = *tmp_p;

        return complex_type(tmp.x,tmp.y);
    }

    CUDA
    complex_type get_lin_SoA(index_type index, size_t box_id)
    {

        const complex_type* ptr = array[index];
        const REAL2* tmp_p = reinterpret_cast<const REAL2*>(ptr + box_id);
        const REAL2 tmp = (*tmp_p);

        return complex_type(tmp.x,tmp.y);
    }

    CUDA
    complex_type* get_lin_SoA_ptr(size_t index, size_t box_id)
    {

        complex_type* ptr = array[index];
        complex_type* val = &ptr[box_id];

        return val;
    }

    CUDA
    complex_type* get_SoA_ptr(index_type l, index_type m, size_t box_id)
    {

        complex_type* ptr = array[offset(l, m)];
        complex_type* val = &ptr[box_id];

        return val;
    }

    CUDA
    complex_type& get_lin_SoA_val(size_t index, size_t box_id)
    {
        complex_type* ptr = array[index];
        return ptr[box_id];
    }

    CUDA
    complex_type getSoA_full(index_type l, index_type m, size_t box_id)
    {
        index_type abs_m = (m<0) ? -m:m;
        const complex_type* ptr = array[offset(l, abs_m)];
        const REAL2* tmp_p = reinterpret_cast<const REAL2*>(ptr + box_id);
        const REAL2 tmp = (*tmp_p);
        complex_type val(tmp.x,tmp.y);
        return (m<0) ? cuda_toggle_sign_if_odd(abs_m, conj(val)) : val;
    }

    CUDA
    value_type get_vectorized_upper(index_type l, index_type m) const
    {
        const REAL2* tmp_p = reinterpret_cast<const REAL2*>(array + offset(l, m));
        const REAL2 tmp = *tmp_p;
        return value_type(tmp.x,tmp.y);
    }
    CUDA
    value_type get_vectorized_lower(index_type l, index_type m) const
    {
        const REAL2* tmp_p = reinterpret_cast<const REAL2*>(array + offset(l, -m));
        const REAL2 tmp = *tmp_p;
        value_type ret(tmp.x,tmp.y);
        return cuda_toggle_sign_if_odd(-m, conj(ret));
    }

    DEVICE
    value_type get_vectorized(size_t index) const
    {
        const REAL2* tmp_p = reinterpret_cast<const REAL2*>(array + index);
        const REAL2 tmp = *tmp_p;
        return value_type(tmp.x,tmp.y);
    }

    DEVICE
    value_type get_vectorized(index_type l, index_type m) const
    {
        return (m < 0) ? get_vectorized_lower(l, m) : get_vectorized_upper(l, m);
    }

    // 0 <= l <= p, 0 <= m <= l
    CUDA
    const value_type & operator () (index_type l, index_type m) const
    {
        return get_upper(l, m);
    }

    // 0 <= l <= p, 0 <= m <= l
    CUDA
    value_type & operator () (index_type l, index_type m)
    {
        return get_upper(l, m);
    }

    CUDA
    value_type & operator () (size_t index)
    {
        return array[index];
    }

    CUDA
    const value_type & operator () (size_t index) const
    {
        return array[index];
    }

    CUDA
    void populate_lower()
    { }  // no-op

    CUDA
    void realloc(value_type * arr){
        array = arr;
    }

    CUDA
    UpperTriangularArray & operator *= ( value_type lambda )
    {
        for (size_type i = 0; i < size(); ++i)
            array[i] *= lambda;

        return *this;
    }

    CUDA
    UpperTriangularArray & operator += ( UpperTriangularArray & other )
    {
        for (size_type i = 0; i < size(); ++i)
            array[i] += other.array[i];

        return *this;
    }

    CUDA
    UpperTriangularArray & operator -= ( UpperTriangularArray & other )
    {
        for (size_type i = 0; i < size(); ++i)
            array[i] -= other.array[i];

        return *this;
    }

    CUDA
    UpperTriangularArray & operator = ( UpperTriangularArray & other )
    {
        for (size_type i = 0; i < size(); ++i)
            array[i] = other.array[i];

        return *this;
    }


    template <typename UpperTriangularArrayD>
    CUDA
    UpperTriangularArray & recast( UpperTriangularArrayD & other )
    {
        for (size_type i = 0; i < size(); ++i)
        {
            complex_type val = other(i);
            array[i].real((Realtype)val.real());
            array[i].imag((Realtype)val.imag());
        }


        return *this;
    }



private:
    value_type *array;
    allocator_type allocator;
    index_type p_;

    CUDA
    size_type size() const
    {
        return size(p_);
    }

    CUDA
    size_type size(index_type p) const
    {
        // sum(1..p+1): there are (p+1) elements in the longest row/col
        return ((p + 1) * (p + 2)) / 2;
    }

    // 0 <= l <= p, 0 <= m <= l
    CUDA
    index_type offset(index_type l, index_type m) const
    {
        //printf("ORIG l p m  %d %d %d\n",(int)l,(int)p_,(int)m);
        //TRIANGULAR_ARRAY_RANGE_ASSERT(0 <= l && l <= p_ && 0 <= m && m <= l);
        //TRIANGULAR_ARRAY_RANGE_ASSERT(0 <= l);
        //TRIANGULAR_ARRAY_RANGE_ASSERT(l <= p_);
        //TRIANGULAR_ARRAY_RANGE_ASSERT(0 <= m);
        //TRIANGULAR_ARRAY_RANGE_ASSERT(m <= l);
        return (l * (l + 1)) / 2 + m;
    }
};



template <typename CoefficientMatrix>
void dump_upper_triangle(const CoefficientMatrix & a, size_t p)
{
    std::cout << std::scientific;
    //std::cout << std::scientific << setprecision(12) << setw() <<
    for (size_t m = p; m <= p; --m) {
        for (size_t l = 0; l <= p; ++l) {
            if (m <= l)
                std::cout << std::setprecision(10) << std::setw(21) << a(l, m) << " ";
            //else
            //std::cout << std::setw(21) << "() ";
        }
        std::cout << std::endl;
    }
}


template <typename CoefficientMatrix>

void dump_upper_triangle_list(const CoefficientMatrix & a, size_t p)
{
    std::cout << "# l\tm\tre\tim" << std::endl;
    for (size_t l = 0; l <= p; ++l) {
        for (size_t m = 0; m <= l; ++m) {
            std::cout << l << "\t" << m << "\t" << a(l, m).real() << "\t" << a(l, m).imag() << std::endl;
        }
    }
}

template <typename CoefficientMatrix>
void dump_all(const CoefficientMatrix & a, size_t p)
{
    std::cout << std::scientific;
    //std::cout << std::scientific << setprecision(12) << setw() <<
    for (ssize_t l = 0; l <= p; ++l){
    for (ssize_t m = -l; m <= l; ++m)
    {
            std::cout << std::setprecision(10) << std::setw(21) << a.get(l, m) << " ";
        //else
        //std::cout << std::setw(21) << "() ";
    }
    std::cout << std::endl;
    }
}

template <typename CoefficientMatrix>
void dump(const CoefficientMatrix & a, size_t p)
{
    std::cout << std::scientific << std::setprecision(10);

    for (ssize_t l = 0; l <= p; ++l)
    {
        for (ssize_t m = 0; m <= l; ++m)
        {
            std::cout << a.get(l, m) << " ";
        }
        std::cout << std::endl;
    }
}

template <typename CoefficientMatrix>
void dump_all_withsignbits( CoefficientMatrix & a, size_t p)
{
    std::cout << std::scientific;
    unsigned int index = 0;
    for (ssize_t l = 0; l <= p; ++l)
    {
        for (ssize_t m = -l; m <= l; ++m)
        {
                std::cout << std::setprecision(10) << std::setw(21) << a.get(l, m) << " ";
                printf("(%lu ",a.bitset.get(index));
                printf("%lu)",a.bitset.get(++index));
                index++;
        }
        std::cout << std::endl;
    }
}

template <typename T, typename architecture>
class UpperTriangularArray_Memory : public Managed<architecture>
{

public:
    typedef T value_type;
    typedef typename architecture::template rebind<value_type>::other allocator_type;
    typedef size_t size_type;
    typedef ssize_t index_type;
    
    UpperTriangularArray_Memory(){}
    
    UpperTriangularArray_Memory(index_type p, size_t num_of_Arrays) : p_(p),vec_size(num_of_Arrays)
    {
        architecture::custom_alloc(array, vec_size*sizeof(value_type*));
        pointer = allocator.allocate(size()*vec_size);
        for (index_type i=0;i<vec_size;i++)
        {
            array[i] = pointer + i*size();
            for (size_type j = 0; j < size(); ++j)
            {
                allocator.construct(&(array[i][j]), value_type());
            }
        }
    }
    
    ~UpperTriangularArray_Memory()
    {
        for (index_type i=0;i<vec_size;i++)
        {
            array[i] = pointer + i*size();
            for (size_type j = 0; j < size(); ++j)
            {
                allocator.destroy(&(array[i][j]));
            }
        }
        allocator.deallocate(pointer,size()*vec_size);
        architecture::custom_free(array);
    }

    //mapping memory from Vector_of_UpperTriangularArray to UpperTriangularArray
    CUDA
    value_type *get_raw_pointer(index_type index){
        return array[index];
    }
    
private:
    value_type* pointer;
    value_type ** array;
    allocator_type allocator;
    index_type p_;
    size_t vec_size;

    CUDA
    size_type size() const
    {
        return size(p_);
    }

    CUDA
    size_type size(index_type p) const
    {
        return ((p + 1) * (p + 2)) / 2;
    }

    CUDA
    index_type offset(index_type l, index_type m) const
    {
        return (l * (l + 1)) / 2 + m;
    }
};



















































// triangular matrix with corners (0,0), (-p,-p), (p,p)
// optimized for row-major access
// only the elements in the upper triangle (0,0), (p,0), (p,p) are actually
// stored and efficiently accessible
template <typename T, typename ComputeLower, typename Allocator = std::allocator<T> >
class UpperTriangularArrayRowMajor
{
    //UpperTriangularArrayRowMajor(const UpperTriangularArrayRowMajor &) = delete;
    //UpperTriangularArrayRowMajor & operator=(const UpperTriangularArrayRowMajor &) = delete;
public:
    typedef T value_type;
    typedef ComputeLower compute_lower_functor;
    typedef typename Allocator::template rebind<value_type>::other allocator_type;
    typedef size_t size_type;
    typedef ssize_t index_type;

    typedef T* row_iterator;

    UpperTriangularArrayRowMajor(index_type p) : p_(p)
    {
        array = allocator.allocate(size());
        for (size_type i = 0; i < size(); ++i)
            allocator.construct(&array[i], value_type());
    }

    ~UpperTriangularArrayRowMajor()
    {
        for (size_type i = 0; i < size(); ++i)
            allocator.destroy(&array[i]);
        allocator.deallocate(array, size());
    }
    CUDA
    size_type p() const
    {
        return p_;
    }
    CUDA
    void fill(const value_type & value = value_type())
    {
        std::fill(array, array + size(), value);
    }
    CUDA
    void zero()
    {
        fill();
    }
    CUDA
    void zero(index_type p)
    {
        std::fill(array, array + size(p), value_type());
    }

    // 0 <= l <= p, 0 <= m <= l
    CUDA
    const value_type & get_upper(index_type l, index_type m) const
    {
        return array[offset(l, m)];
    }

    // 0 <= l <= p, 0 <= m <= l
    CUDA
    value_type & get_upper(index_type l, index_type m)
    {
        return array[offset(l, m)];
    }

    // 0 <= l <= p, -l <= m <= -1
    CUDA
    value_type get_lower(index_type l, index_type m) const
    {
        return compute_lower_functor()(l, m, *this);
    }

    // 0 <= l <= p, -l <= m <= l
    CUDA
    value_type get(index_type l, index_type m) const
    {
        return (m < 0) ? get_lower(l, m) : get_upper(l, m);
    }

    // 0 <= l <= p, 0 <= m <= l
    CUDA
    const value_type & operator () (index_type l, index_type m) const
    {
        return get_upper(l, m);
    }

    // 0 <= l <= p, 0 <= m <= l
    CUDA
    value_type & operator () (index_type l, index_type m)
    {
        return get_upper(l, m);
    }
    CUDA
    row_iterator get_row_iterator(index_type l, index_type m)
    {
        return array + offset(l, m);
    }

    void populate_lower()
    { }  // no-op

private:
    value_type * array;
    allocator_type allocator;
    index_type p_;

    CUDA
    size_type size() const
    {
        return size(p_);
    }
    CUDA
    size_type size(index_type p) const
    {
        // sum(1..p+1): there are (p+1) elements in the longest row/col
        return ((p + 1) * (p + 2)) / 2;
    }

    // 0 <= l <= p, 0 <= m <= l
    CUDA
    index_type offset(index_type l, index_type m) const
    {
        //TRIANGULAR_ARRAY_RANGE_ASSERT(0 <= l && l <= p_ && 0 <= m && m <= l);
        return m * (2 * p_ - m + 1) / 2 + l;
    }
};


// triangular matrix with corners (0,0), (-p,-p), (p,p)
// for elements at positions (l,m) with 0 <= l <= p, -l <= m <= l
// optimized for column-major access
template <typename T, typename ComputeLower, typename Allocator = std::allocator<T> >
class TriangularArray
{
    //TriangularArray(const TriangularArray &) = delete;
    //TriangularArray & operator=(const TriangularArray &) = delete;
public:
    typedef T value_type;
    typedef ComputeLower compute_lower_functor;
    typedef size_t size_type;

    typedef typename Allocator::template rebind<value_type>::other allocator_type;

    typedef ssize_t index_type;
    typedef const T * const_column_iterator;

    TriangularArray(index_type p) : p_(p)
    {
        array = allocator.allocate(size());
        for (size_type i = 0; i < size(); ++i)
            allocator.construct(&array[i], value_type());

        //array = new value_type[size()]();
    }

    ~TriangularArray()
    {
        for (size_type i = 0; i < size(); ++i)
            allocator.destroy(&array[i]);

        allocator.deallocate(array, size());
        //delete[] array;
    }

    CUDA
    size_type p() const
    {
        return p_;
    }
    CUDA
    void fill(const value_type & value = value_type())
    {
        std::fill(array, array + size(), value);
    }
    CUDA
    void zero()
    {
        fill();
    }
    CUDA
    void zero(index_type p)
    {
        std::fill(array, array + size(p), value_type());
    }

    // 0 <= l <= p, 0 <= m <= l
    CUDA
    const value_type & get_upper(index_type l, index_type m) const
    {
        return get(l, m);
    }

    // 0 <= l <= p, 0 <= m <= l
    CUDA
    value_type & get_upper(index_type l, index_type m)
    {
        return get(l, m);
    }

    // 0 <= l <= p, -l <= m <= -1
    CUDA
    const value_type & get_lower(index_type l, index_type m) const
    {
        return get(l, m);
    }

    // 0 <= l <= p, -l <= m <= -1
    CUDA
    value_type & get_lower(index_type l, index_type m)
    {
        return get(l, m);
    }

    // 0 <= l <= p, -l <= m <= l
    CUDA
    const value_type & get(index_type l, index_type m) const
    {
        return array[offset(l, m)];
    }

    // 0 <= l <= p, -l <= m <= l
    CUDA
    value_type & get(index_type l, index_type m)
    {
        return array[offset(l, m)];
    }

    // 0 <= l <= p, 0 <= m <= l
    CUDA
    const value_type & operator () (index_type l, index_type m) const
    {

        return get_upper(l, m);
    }

    // 0 <= l <= p, 0 <= m <= l
    CUDA
    value_type & operator () (index_type l, index_type m)
    {
        return get_upper(l, m);
    }
    CUDA
    void populate_lower()
    {
        compute_lower_functor compute_lower;
        for (index_type l = 1; l <= p_; ++l)
            for (index_type m = -1; m >= -l; --m)
                array[offset(l, m)] = compute_lower(l, m, *this);
    }
    CUDA
    const_column_iterator get_column_iterator(index_type l, index_type m) const
    {
        return array + offset(l, m);
    }

private:
    value_type * array;
    index_type p_;
    allocator_type allocator;

    CUDA
    size_type size() const
    {
        return size(p_);
    }
    CUDA
    size_type size(index_type p) const
    {
        // upper (0,0), (p,0), (p,p): sum(1..p+1)
        // lower (-1,-1), (-p,-p), (p,-1): sum(1..p)
        return (p + 1) * (p + 1);
    }

    // 0 <= l <= p, 0 <= m <= l
    CUDA
    index_type offset(index_type l, index_type m) const
    {
        //TRIANGULAR_ARRAY_RANGE_ASSERT(0 <= l && l <= p_ && -l <= m && m <= l);
        return l * l + l + m;
    }
};

}//namespace end

#endif
// vim: et:ts=4:sw=4
