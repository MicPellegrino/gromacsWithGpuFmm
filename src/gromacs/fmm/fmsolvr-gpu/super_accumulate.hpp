#ifndef _BK_FMM_super_accumulate_hpp_
#define _BK_FMM_super_accumulate_hpp_

#include <functional>

namespace gmx_gpu_fmm{

template <typename Storage, typename SuSI>
class super_accumulate
{
    typedef Storage storage_type;
    typedef typename storage_type::value_type storage_value_type;
public:
    typedef SuSI value_type;
    typedef typename storage_type::index_type index_type;

    super_accumulate(storage_type & storage) : storage(storage)
    { }

    void operator () (index_type l, index_type m, const value_type & value)
    {
        storage(l, m) += eindampfen(value);
    }

private:
    template <typename T>
    storage_value_type eindampfen(const typename COMPLEX_GENERATOR<T>::type & value)
    {
        Vec_traits<T> VT;
        return storage_value_type(VT.supersum(value.real()), VT.supersum(value.imag()));
    }

    storage_value_type NIXeindampfen(const value_type & value)
    {
        Vec_traits<value_type> VT;
        return VT.supersum(value);
    }

    storage_type & storage;
};

template <typename StorageT, typename SuSI>
class just_complex_super_accumulate : public std::binary_function<StorageT, SuSI, void>
{
    typedef StorageT storage_type;
public:
    typedef SuSI value_type;

    void operator () (storage_type & target, const value_type & value) const
    {
        Vec_traits<typename value_type::value_type> VT;
        target += storage_type(VT.supersum(value.real()), VT.supersum(value.imag()));;
    }
};

template <typename StorageT, typename SuSI>
class just_complex_interleaved_super_accumulate : public std::binary_function<StorageT, SuSI, void>
{
    typedef StorageT storage_type;
public:
    typedef SuSI value_type;

    void operator () (storage_type & target, const value_type & value) const
    {
        typedef typename value_type::value_type super_simd_type;
        //typedef typename storage_type::value_type simd_type;
        typedef storage_type simd_type;
        Vec_traits<super_simd_type> SuVT;
        Vec_traits<simd_type> VT;
        simd_type re = SuVT.supersum(value.real());
        simd_type im = SuVT.supersum(value.imag());
        simd_type interleaved = VT.interleavedsumpair(re, im);
        //simd_type re_sult = target.real() + interleaved;
        //target.real(re_sult);
        target += interleaved;
    }
};

template <typename Storage, typename SuSI>
class scalar_accumulate
{
    typedef Storage storage_type;
    typedef typename storage_type::value_type storage_value_type;
public:
    typedef SuSI value_type;
    typedef typename storage_type::index_type index_type;

    scalar_accumulate(storage_type & storage) : storage(storage)
    { }

    void operator () (index_type l, index_type m, const value_type & value)
    {
#if 1
        storage(l, m) += eindampfen(value);
#else
        storage_value_type eingedampft = eindampfen(value);
        __m128d val = *reinterpret_cast<__m128d *>(&eingedampft);
        __m128d * res = reinterpret_cast<__m128d *>(&storage(l, m));
        *res += val;
#endif
    }

private:
    template <typename T>
    storage_value_type eindampfen(const typename COMPLEX_GENERATOR<T>::type & value)
    {
        Vec_traits<T> SuVT;
        typedef typename Vec_traits<T>::nested_simd_type nested_simd_type;
        Vec_traits<nested_simd_type> VT;
        //return storage_value_type(VT.sum(SuVT.supersum(value.real())), VT.sum(SuVT.supersum(value.imag())));
        nested_simd_type re = SuVT.supersum(value.real());
        nested_simd_type im = SuVT.supersum(value.imag());
        nested_simd_type interleaved_sum = VT.sumpair(re, im);
#if 1
        return storage_value_type(VT.element(interleaved_sum, 0), VT.element(interleaved_sum, 1));
#else
        storage_value_type res;
        *reinterpret_cast<__m128d *>(&res) = _mm256_castpd256_pd128(interleaved_sum);
        return res;
#endif
    }

    storage_type & storage;
};

template <typename StorageT, typename SuSI>
class just_complex_scalar_accumulate : public std::binary_function<StorageT, SuSI, void>
{
    typedef StorageT storage_type;
public:
    typedef SuSI value_type;

    void operator () (storage_type & target, const value_type & value) const
    {
#if 1
        target += eindampfen(value);
#else
        storage_type eingedampft = eindampfen(value);
        __m128d val = *reinterpret_cast<__m128d *>(&eingedampft);
        __m128d * res = reinterpret_cast<__m128d *>(&target);
        *res += val;
#endif
    }

private:
    template <typename T>
    storage_type eindampfen(const typename COMPLEX_GENERATOR<T>::type & value) const
    {
        Vec_traits<T> SuVT;
        typedef typename Vec_traits<T>::nested_simd_type nested_simd_type;
        Vec_traits<nested_simd_type> VT;
        //return storage_value_type(VT.sum(SuVT.supersum(value.real())), VT.sum(SuVT.supersum(value.imag())));
        nested_simd_type re = SuVT.supersum(value.real());
        nested_simd_type im = SuVT.supersum(value.imag());
        nested_simd_type interleaved_sum = VT.sumpair(re, im);
#if 1
        return storage_type(VT.element(interleaved_sum, 0), VT.element(interleaved_sum, 1));
#else
        storage_type res;
        *reinterpret_cast<__m128d *>(&res) = _mm256_castpd256_pd128(interleaved_sum);
        return res;
#endif
    }

};

}//namespace end

#endif
