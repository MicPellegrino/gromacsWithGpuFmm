#ifndef _BK_FMM_parity_sign_hpp_
#define _BK_FMM_parity_sign_hpp_

namespace gmx_gpu_fmm{

template <typename Real>
CUDA
__forceinline__
Real powminus1(size_t x)
{
    return (x & 1) ? Real(-1.) : Real(1);
}

template <typename Integral>
CUDA
__forceinline__
bool is_odd(Integral x)
{
    return x & 1;
}

// (-1)^k * val
template <typename T, typename Integral>
CUDA
__forceinline__
T toggle_sign_if_odd(Integral k, T val)
{
    return is_odd(k) ? -val : val;
}

template <typename TT, typename T, typename Integral>
CUDA
__forceinline__
T toggle_sign_if_odd_realpart(Integral k, T val)
{
    T realval = s_odd(k) ? -val.real() : val.real();
    T imagval = val.imag();
    TT ret(realval,imagval);
    return ret;
}

template <typename Integral>
DEVICE
__forceinline__
#ifndef GMX_FMM_DOUBLE
    int bitmask(Integral x)
#else
    long long int bitmask(Integral x)
#endif

{
#ifndef GMX_FMM_DOUBLE
    return (x & 1)<<31;
#else
    return (x & 1)<<63;
#endif
}

template <typename T, typename Integral>
DEVICE
__forceinline__
T cuda_toggle_sign_if_odd(Integral k, T val)
{

#ifdef __CUDACC__

    double x = val.real();
    double y = val.imag();

#ifndef GMX_FMM_DOUBLE
    int mask = bitmask(k);
    T v(__int_as_float( __float_as_int(x)^mask),__int_as_float( __float_as_int(y)^mask));
#else
    long long int mask = bitmask(k);
    T v(__longlong_as_double( __double_as_longlong(x)^mask),__longlong_as_double( __double_as_longlong(y)^mask));
#endif
    return v;
#else

    return toggle_sign_if_odd(k, val);

#endif

}

template <typename T>
DEVICE
__forceinline__
T cuda_toggle_sign_mask(long long int mask, T &val)
{

#ifdef __CUDACC__
    double x = val.real();
    double y = val.imag();
    T v(__longlong_as_double( __double_as_longlong(x)^mask),__longlong_as_double( __double_as_longlong(y)^mask));
    return v;
#endif

}

template <typename T>
DEVICE
__forceinline__
T cuda_toggle_sign_mask(int mask, T &val)
{

#ifdef __CUDACC__
    float x = val.real();
    float y = val.imag();
    T v(__int_as_float( __float_as_int(x)^mask),__int_as_float( __float_as_int(y)^mask));
    return v;
#endif

}

}//namespace end


#endif
