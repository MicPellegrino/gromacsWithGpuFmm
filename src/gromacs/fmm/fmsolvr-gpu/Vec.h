#ifndef _BK_FMM_Vec_h_
#define _BK_FMM_Vec_h_

#include "Vec_base.h"
#include "Vec_Array.h"
#include "Vec_super.h"

namespace gmx_gpu_fmm
{


//#ifdef __NVCC__
//#undef __AVX__
//#endif

template <typename Real>
struct simd_engine_none
{
	typedef Real type;
};

template <typename Real>
struct simd_engine_array8
{
	typedef TrivialArray<Real, 8> type;
};

#ifdef __AVX__
template <typename Real>
struct simd_engine_avx
{

};

template <>
struct simd_engine_avx<float>
{
    //typedef __m256 type;
    typedef float type;
};

template <>
struct simd_engine_avx<double>
{
    //typedef __m256d type;
    typedef double type;
};
#endif  // __AVX__

#ifdef __SSE__
template <typename Real>
struct simd_engine_sse
{

};

template <>
struct simd_engine_sse<float>
{
    //typedef __m128 type;
    typedef float type;
};
#endif  // __SSE__

#ifdef __SSE2__
template <>
struct simd_engine_sse<double>
{
    //typedef __m128d type;
    typedef double type;
};
#endif  // __SSE2__
//#endif

template <typename Real>
struct simd_engine_default
{

#if defined (__AVX__)
	typedef typename simd_engine_avx<Real>::type type;
#elif defined (__SSE__)
	typedef typename simd_engine_sse<Real>::type type;
#else
	typedef typename simd_engine_none<Real>::type type;
#endif
};

template <>
struct simd_engine_default<long double>
{
	typedef typename simd_engine_none<long double>::type type;
};


//#if defined(_GLIBCXX_USE_FLOAT128)
//template <>
//struct simd_engine_default<__float128>
//{
//	typedef typename simd_engine_none<__float128>::type type;
//};
//#endif

}//namespace end
#endif
