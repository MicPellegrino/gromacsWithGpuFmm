#ifndef XYZQ_HPP
#define XYZQ_HPP

#include <iostream>
#include "cuda_keywords.hpp"

namespace gmx_gpu_fmm{

template <typename Real>
class XYZQ {

public:

    typedef Real value_type;
    value_type x, y, z, q;

    CUDA
    XYZQ();

    CUDA
    XYZQ(value_type x, value_type y, value_type z, value_type q);

};

template <typename Real>
CUDA
XYZQ<Real> make_xyzq(const Real & x, const Real & y, const Real & z, const Real & q)
{
    return XYZQ<Real>(x, y, z, q);
}

template <typename Real, typename XYZ>
CUDA
XYZQ<Real> make_xyzq(const XYZ & xyz, const Real & q)
{
    return XYZQ<Real>(xyz.x, xyz.y, xyz.z, q);
}

template <typename Real>
std::ostream & operator << (std::ostream & os, const XYZQ<Real> & xyzq)
{
    os << "(" << xyzq.x << ", " << xyzq.y << ", " << xyzq.z << "; " << xyzq.q<< ")";
    return os;
}

}//namespace end

#endif // XYZQ_HPP
