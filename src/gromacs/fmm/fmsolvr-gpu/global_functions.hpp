#ifndef GLOBAL_FUNCTIONS_HPP
#define GLOBAL_FUNCTIONS_HPP

#include <math.h>
#include <vector>
#include <algorithm>
#include <cassert>
#include <stdlib.h>
#include <fstream>
#include "linear_algebra.hpp"
#include "xyzq.hpp"
#include "xyz.hpp"

namespace gmx_gpu_fmm{

extern int get_env_int(int dflt, const char* vrnm);

extern size_t boxes_above_depth(size_t d);

extern size_t boxes_on_depth(size_t d);

extern size_t make_boxid(size_t x, size_t y, size_t z, unsigned depth);

template <typename Real>
inline Real reciprocal(Real a)
{
    return  Real(1.) / a;
}

template <typename Real>
inline Real rsqrt(Real r)
{
    return Real(1.) / std::sqrt(r);
}

template <typename Real33>
Real33 change_of_basis_from_standard_basis(const Real33 & abc)
{
    typedef typename Real33::value_type Real3;
    typedef typename Real3::value_type Real;


    /*
    linear_algebra::SqMatrix<Real> B(3), Binv(3);
    B[0][0] = abc.a.x;
    B[0][1] = abc.a.y;
    B[0][2] = abc.a.z;
    B[1][0] = abc.b.x;
    B[1][1] = abc.b.y;
    B[1][2] = abc.b.z;
    B[2][0] = abc.c.x;
    B[2][1] = abc.c.y;
    B[2][2] = abc.c.z;
    if (!linear_algebra::inverse(3, B, Binv))
        throw "basis matrix not invertible - not a valid basis";

    Real33 result = Real33(Real3(Binv[0][0], Binv[0][1], Binv[0][2]),
                           Real3(Binv[1][0], Binv[1][1], Binv[1][2]),
                           Real3(Binv[2][0], Binv[2][1], Binv[2][2]));

    std::cout<<result.a<<std::endl;
    */

    //HACK HACK NUMERICS ISSUES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //removes the problem of particles placed at the edge of the simulation box
    //the particles with positon x == simulation_box_size are put in the box as well
    Real size = reciprocal(std::nextafter(abc.a.x, abc.a.x + 1));
    //Real size = reciprocal(abc.a.x);
    //printf("size %.30e %.30e\n", 1.0/abc.a.x,size);
    //Real size = abc.a.x;
    Real33 result = Real33(Real3(size, 0., 0.),
                           Real3(0., size, 0.),
                           Real3(0., 0., size));


    return result;
}


template <typename Real33>
Real33 change_of_basis_to_standard_basis(const Real33 & abc)
{
    //std::cout<<abc.a<<abc.b<<abc.c<<std::endl;
    return abc;
}

template <typename RealOut, typename RealIn = RealOut>
struct fake_AoS
{
    typedef RealOut Real;
    typedef const RealIn* const_pointer;
    typedef XYZQ<Real> Real4;
    typedef Real4 value_type;
    typedef RealIn vec[3];

    const_pointer x_, y_, z_, q_;
    size_t n;

    fake_AoS(const_pointer x, const_pointer y, const_pointer z, const_pointer qq, size_t n)
        : x_(x), y_(y), z_(z), q_(qq), n(n)
    { }

    Real4 operator () (size_t i) const
    {
        return make_xyzq<Real>(x_[i], y_[i], z_[i], q_[i]);
    }

    const RealOut& x(size_t i) const
    {
        return x_[i];
    }

    const RealOut& y(size_t i) const
    {
        return y_[i];
    }

    const RealOut& z(size_t i) const
    {
        return z_[i];
    }

    const RealOut& q(size_t i) const
    {
        return q_[i];
    }
};

template <typename RealOut, typename RealIn = RealOut>
struct fake_AoS_GMX
{
    typedef RealOut Real;
    typedef const RealIn* const_pointer;
    typedef XYZQ<Real> Real4;
    typedef Real4 value_type;
    typedef RealIn vec[3];

    const vec* xyz;
    const_pointer q_;
    size_t n;

    fake_AoS_GMX(const vec* xyz, const RealIn* qq, const size_t n)
        : xyz(xyz), q_(qq), n(n)
    { }

    Real4 operator () (size_t i) const
    {
        return make_xyzq<Real>(xyz[i][0], xyz[i][1], xyz[i][2], q_[i]);
    }


    template <typename Real4Compatible>
    void set(Real4Compatible& v, Real scale) const
    {
        for (int i = 0; i < n; ++i)
        {
            v[i].x = xyz[i][0] * scale;
            //assert(v[i].x >= 0);
            v[i].y = xyz[i][1] * scale;
            //assert(v[i].y >= 0);
            v[i].z = xyz[i][2] * scale;
            //assert(v[i].z >= 0);
            v[i].q = q_[i];
        }
    }

    const RealOut x(size_t i) const
    {
        return xyz[i][0];
    }

    const RealOut y(size_t i) const
    {
        return xyz[i][1];
    }

    const RealOut z(size_t i) const
    {
        return xyz[i][2];
    }

    const RealOut q(size_t i) const
    {
        return q_[i];
    }
};

}//namespace end

#endif // GLOBAL_FUNCTIONS_HPP
