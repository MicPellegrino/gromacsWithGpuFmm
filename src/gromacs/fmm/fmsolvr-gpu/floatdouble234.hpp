#ifndef FLOATDOUBLE234_HPP
#define FLOATDOUBLE234_HPP

#include "cuda_keywords.hpp"
#include "data_type.hpp"

namespace gmx_gpu_fmm{

CUDA __forceinline__
REAL2 & operator +=(REAL2& a, const REAL2& b)
{
    a.x += b.x;
    a.y += b.y;
    return a;
}

CUDA __forceinline__
REAL3 & operator +=(REAL3& a, const REAL3& b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

CUDA __forceinline__
REAL4 & operator +=(REAL4& a, const REAL4& b)
{
    a.x +=b.x;
    a.y +=b.y;
    a.z +=b.z;
    a.w +=b.w;
    return a;
}

CUDA __forceinline__
REAL4 & operator +=(REAL4& a, const REAL3& b)
{
    a.x +=b.x;
    a.y +=b.y;
    a.z +=b.z;
    return a;
}

CUDA __forceinline__
REAL2 & operator -=(REAL2& a, const REAL2& b)
{
    a.x -= b.x;
    a.y -= b.y;
    return a;
}

CUDA __forceinline__
REAL3 & operator -=(REAL3& a, const REAL3& b)
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    return a;
}

CUDA __forceinline__
REAL4 & operator -=(REAL4& a, const REAL4& b)
{
    a.x -=b.x;
    a.y -=b.y;
    a.z -=b.z;
    a.w -=b.w;
    return a;
}

CUDA __forceinline__
REAL4 & operator -=(REAL4& a, const REAL3& b)
{
    a.x -=b.x;
    a.y -=b.y;
    a.z -=b.z;
    return a;
}

CUDA __forceinline__
REAL2  operator + (REAL2 a, const REAL2& b)
{
    a +=b;
    return a;
}

CUDA __forceinline__
REAL3  operator + (REAL3 a, const REAL3& b)
{
    a +=b;
    return a;
}

CUDA __forceinline__
REAL4  operator + (REAL4 a, const REAL4& b)
{
    a +=b;
    return a;
}


CUDA __forceinline__
REAL2  operator - (REAL2 a, const REAL2& b)
{
    a -=b;
    return a;
}

CUDA __forceinline__
REAL3  operator - (REAL3 a, const REAL3& b)
{
    a -=b;
    return a;
}

CUDA __forceinline__
REAL4  operator - (REAL4 a, const REAL4& b)
{
    a -=b;
    return a;
}

CUDA __forceinline__
REAL2 & operator *=(REAL2& a, const REAL2& b)
{
    a.x *= b.x;
    a.y *= b.y;
    return a;
}

CUDA __forceinline__
REAL3 & operator *=(REAL3& a, const REAL3& b)
{
    a.x *= b.x;
    a.y *= b.y;
    a.z *= b.z;
    return a;
}

CUDA __forceinline__
REAL4 & operator *=(REAL4& a, const REAL4& b)
{
    a.x *=b.x;
    a.y *=b.y;
    a.z *=b.z;
    a.w *=b.w;
    return a;
}

CUDA __forceinline__
REAL2  operator * (REAL2 a, const REAL2& b)
{
    a *=b;
    return a;
}

CUDA __forceinline__
REAL3  operator * (REAL3 a, const REAL3& b)
{
    a *=b;
    return a;
}

CUDA __forceinline__
REAL4  operator * (REAL4 a, const REAL4& b)
{
    a *=b;
    return a;
}

CUDA __forceinline__
REAL2 & operator *=(REAL2& a, const REAL& b)
{
    a.x *= b;
    a.y *= b;
    return a;
}

CUDA __forceinline__
REAL3 & operator *=(REAL3& a, const REAL& b)
{
    a.x *= b;
    a.y *= b;
    a.z *= b;
    return a;
}

CUDA __forceinline__
REAL4 & operator *=(REAL4& a, const REAL& b)
{
    a.x *=b;
    a.y *=b;
    a.z *=b;
    a.w *=b;
    return a;
}

CUDA __forceinline__
REAL2  operator * (REAL2 a, const REAL& b)
{
    a *=b;
    return a;
}

CUDA __forceinline__
REAL3  operator * (REAL3 a, const REAL& b)
{
    a *=b;
    return a;
}

CUDA __forceinline__
REAL4  operator * (REAL4 a, const REAL& b)
{
    a *=b;
    return a;
}

}//namespace end

#endif
