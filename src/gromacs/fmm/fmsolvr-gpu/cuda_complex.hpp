#ifndef CUDA_COMPLEX_HPP
#define CUDA_COMPLEX_HPP

#include "cuda_atomics.hpp"
#include "data_type.hpp"
#include "cuda_keywords.hpp"
#include <math.h>
#include <sstream>

namespace gmx_gpu_fmm{

/*! \brief Implements the complex number functionality for CUDA.
 *  \tparam T Scalar data type.
 */
template<class T>
class  complex
{

public:
    //! Underlying value type
    typedef T value_type;

protected:
    //! Real part
    value_type __re_;
    //! Imaginary part
    value_type __im_;
public:


    /*!
     * \brief Constructor.
     * \param __re  Initial value for the real part.
     * \param __im  Initial value for the imaginary part.
     */
    CUDA
    complex(const value_type& __re = value_type(), const value_type& __im = value_type())
        :
          __re_(__re), __im_(__im)
    {}

    /*!
     * \brief Copy constructor.
     * \param __c   Complex initial value.
     */
    template<class X>
    CUDA
    complex(const complex<X>& __c)
        :
          __re_(__c.real()), __im_(__c.imag()) {}


    /*!
     * \brief Getter function
     * \return Real part
     */
    CUDA
    value_type real() const
    {
        return __re_;
    }

    /*!
     * \brief Getter function
     * \return Imaginary part part
     */
    CUDA
    value_type imag() const
    {
        return __im_;
    }

    /*!
     * \brief Setter function for real part
     */
    CUDA
    void real(value_type __re)
    {
        __re_ = __re;
    }

    /*!
     * \brief Setter function for imaginary part
     */
    CUDA
    void imag(value_type __im)
    {
        __im_ = __im;
    }

    /*!
     * \brief Implementation of the operator = for real part
     * \param __re  Real type
     */
    CUDA
    complex& operator = (const value_type& __re)
    {
        __re_ = __re;
        __im_ = value_type();
        return *this;
    }

    /*!
     * \brief Implementation of the operator += for real part
     * \param __re  Real type
     */
    CUDA
    complex& operator += (const value_type& __re)
    {
        __re_ += __re;
        return *this;
    }

    CUDA
    complex& operator -= (const value_type& __re)
    {
        __re_ -= __re;
        return *this;
    }

    CUDA
    complex& operator *= (const value_type& __re)
    {
        __re_ *= __re;
        __im_ *= __re;
        return *this;
    }

    CUDA
    complex& operator /= (const value_type& __re)
    {
        __re_ /= __re;
        __im_ /= __re;
        return *this;
    }

    template<class X>
    CUDA complex& operator = (const complex<X>& __c)
    {
        __re_ = __c.real();
        __im_ = __c.imag();
        return *this;
    }

    template<class X>
    CUDA
    complex& operator += (const complex<X>& __c)
    {
        __re_ += __c.real();
        __im_ += __c.imag();
        return *this;
    }

    template<class X>
    CUDA
    complex& operator -= (const complex<X>& __c)
    {
        __re_ -= __c.real();
        __im_ -= __c.imag();
        return *this;
    }

    template<class X>
    CUDA
    complex& operator *= (const complex<X>& __c)
    {
        *this = *this * __c;
        return *this;
    }

    template<class X>
    CUDA
    complex& operator /= (const complex<X>& __c)
    {
        *this = *this / __c;
        return *this;
    }

    template<class X>
    CUDA
    complex& operator () (const value_type& re,const value_type& im)
    {
        __re_ = re;
        __im_ = im;
        return *this;
    }

    template<typename X>
    DEVICE
    __forceinline__
    complex& operator %= (const complex<X>& c)
    {
        const REAL2* tmp_p = reinterpret_cast<const REAL2*>(&c);
        const REAL2 tmp = *tmp_p;

        __atomicAdd(&(__re_), tmp.x);
        __atomicAdd(&(__im_), tmp.y);

        return *this;
    }

    template<typename X>
    DEVICE
    __forceinline__
    void atomic_adder(const complex<X> &c )
    {

        const REAL2* tmp_p = reinterpret_cast<const REAL2*>(&c);
        const REAL2 tmp = *tmp_p;

        __atomicAdd(&(__re_), tmp.x);
        __atomicAdd(&(__im_), tmp.y);
    }
};

template<class T>
CUDA __forceinline__
complex<T> operator + (const complex<T>& __x, const complex<T>& __y)
{
    complex<T> __t(__x);
    __t += __y;
    return __t;
}

template<class T>
CUDA __forceinline__
complex<T> operator - (const complex<T>& __x, const complex<T>& __y)
{
    complex<T> __t(__x);
    __t -= __y;
    return __t;
}

template<class T>
CUDA __forceinline__
complex<T> operator * (const complex<T>& __x, const T& __y)
{
    complex<T> __t(__x);
    __t *= __y;
    return __t;
}

template<class T>
CUDA __forceinline__
complex<T> operator * (const T& __y, const complex<T>& __x)
{
    complex<T> __t(__x);
    __t *= __y;
    return __t;
}

template<class T>
CUDA __forceinline__
complex<T> operator - (const complex<T>& x)
{
    return complex<T>(-x.real(), -x.imag());
}

template<class T>
CUDA __forceinline__
complex<T> operator * (const complex<T>& z, const complex<T>& w)
{

    T a = z.real();
    T b = z.imag();
    T c = w.real();
    T d = w.imag();

    T ac = a * c;
    T bd = b * d;

    T ad = a * d;
    T bc = b * c;

    T x = ac - bd;
    T y = ad + bc;

    return complex<T>(x, y);
}

template<class T>
CUDA __forceinline__
complex<T> conj(const complex<T>& x)
{
    return complex<T>(x.real(), -x.imag());
}

template<class T, class CharT, class Traits>
std::basic_ostream<CharT, Traits>& operator << (std::basic_ostream<CharT, Traits>& os, const complex<T>& x)
{
    os <<"(" << x.real() << ',' << x.imag() <<")";
    return os;
}

}//namespace end

#endif







