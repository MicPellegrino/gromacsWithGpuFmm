#ifndef XYZ_HPP
#define XYZ_HPP

#include <ostream>
#include <cmath>
#include <math.h>
#include "cuda_keywords.hpp"
#include "xyzq.hpp"
#include "data_type.hpp"

namespace gmx_gpu_fmm{

template <typename Real>
class XYZ{

public:

    typedef Real value_type;
    //typedef REAL3 values;
    value_type x, y, z;

    CUDA
    XYZ(value_type xyz);

    CUDA
    XYZ();

    CUDA
    XYZ(value_type x, value_type y, value_type z);

    template <typename XYZCompat>
    CUDA
    explicit XYZ(const XYZCompat& something);

    CUDA
    explicit XYZ(const value_type (& array)[3]);

    CUDA
    explicit XYZ(const XYZQ<value_type> &something);

    CUDA
    XYZ operator -();

    CUDA
    XYZ operator + (const XYZ & a) const;

    CUDA
    XYZ operator - (const XYZ & a) const;

    CUDA
    XYZ operator * (const value_type s) const;

    CUDA
    XYZ operator / (const value_type s) const;

    CUDA
    XYZ & operator += (const XYZ & a);

    CUDA
    XYZ & operator -= (const XYZ & a);

    CUDA
    XYZ & operator *= (const value_type s);

    CUDA
    XYZ & operator *= (const XYZ & a);

    CUDA
    XYZ & operator /= (const XYZ & a);

    CUDA
    XYZ & operator /= (const value_type a);

    CUDA
    XYZ xyz_sqrt() const;

    CUDA
    XYZ xyz_abs() const;

    DEVICE
    XYZ & operator %= (const XYZ & a);
};

template <typename Real>
std::ostream & operator << (std::ostream & os, const XYZ<Real> & xyz)
{
    os << "(" << xyz.x << ", " << xyz.y << ", " << xyz.z << ")";
    return os;
}

}//namespace end

#endif // XYZ_HPP
