#include "xyz.hpp"

namespace gmx_gpu_fmm{

template <typename Real>
CUDA
XYZ<Real>::XYZ(value_type xyz):x(xyz), y(xyz), z(xyz){}

/*
template <typename Real>
CUDA
XYZ<Real>::XYZ(): x(value_type()), y(value_type()), z(value_type())
{ }
*/

template <typename Real>
CUDA
XYZ<Real>::XYZ()
{ }

template <typename Real>
CUDA
XYZ<Real>::XYZ(value_type x, value_type y, value_type z): x(x), y(y), z(z) {}

template <typename Real> template <typename XYZCompat>
CUDA
XYZ<Real>::XYZ(const XYZCompat &something): x(something.x), y(something.y), z(something.z)
{ }

template <typename Real>
CUDA
XYZ<Real>::XYZ(const XYZQ<value_type> &something): x(something.x), y(something.y), z(something.z)
{ }

template <typename Real>
CUDA
XYZ<Real>::XYZ(const value_type (& array)[3]): x(array[0]), y(array[1]), z(array[2])
{ }

template <typename Real>
CUDA
XYZ<Real> XYZ<Real>::operator - ()
{
    return XYZ(-x, -y, -z);
}

template <typename Real>
CUDA
XYZ<Real> XYZ<Real>::operator + (const XYZ & a) const
{
    return XYZ(x + a.x, y + a.y, z + a.z);
}

template <typename Real>
CUDA
XYZ<Real> XYZ<Real>::operator - (const XYZ & a) const
{
    return XYZ(x - a.x, y - a.y, z - a.z);
}

template <typename Real>
CUDA
XYZ<Real> XYZ<Real>::operator * (const value_type s) const
{
    return XYZ(x * s, y * s, z * s);
}

template <typename Real>
CUDA
XYZ<Real> XYZ<Real>::operator / (const value_type s) const
{
    return XYZ(x / s, y / s, z / s);
}

template <typename Real>
CUDA
XYZ<Real>& XYZ<Real>::operator += (const XYZ & a)
{
    x += a.x;
    y += a.y;
    z += a.z;
    return *this;
}

template <typename Real>
CUDA
XYZ<Real>& XYZ<Real>::operator -= (const XYZ & a)
{
    x -= a.x;
    y -= a.y;
    z -= a.z;
    return *this;
}

template <typename Real>
CUDA
XYZ<Real>& XYZ<Real>::operator *= (const value_type s)
{
    x *= s;
    y *= s;
    z *= s;
    return *this;
}

template <typename Real>
CUDA
XYZ<Real>& XYZ<Real>::operator *= (const XYZ & a)
{
    x*=a.x;
    y*=a.y;
    z*=a.z;
    return *this;
}

template <typename Real>
CUDA
XYZ<Real>& XYZ<Real>::operator /= (const XYZ & a)
{
    x/=a.x;
    y/=a.y;
    z/=a.z;
    return *this;
}

template <typename Real>
CUDA
XYZ<Real>& XYZ<Real>::operator /= (const value_type s)
{
    x /= s;
    y /= s;
    z /= s;
    return *this;
}

template <typename Real>
CUDA
XYZ<Real> XYZ<Real>::xyz_sqrt() const
{
    return XYZ( sqrt(x), sqrt(y), sqrt(z) );
}

template <typename Real>
CUDA
XYZ<Real> XYZ<Real>::xyz_abs() const
{
    return XYZ( abs(x), abs(y), abs(z) );
}

template <typename Real>
DEVICE
XYZ<Real>& XYZ<Real>::operator %= (const XYZ & a)
{

#ifdef __CUDACC
    __atomicAdd(&x,a.x);
    __atomicAdd(&y,a.y);
    __atomicAdd(&z,a.z);
#else
    x += a.x;
    y += a.y;
    z += a.z;
#endif

    return *this;
}

template class XYZ<double>;
template class XYZ<float>;


}//namespace end



