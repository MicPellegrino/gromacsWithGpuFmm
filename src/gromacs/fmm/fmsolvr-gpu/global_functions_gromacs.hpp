#ifndef GLOBAL_FUNCTIONS_GROMACS_HPP
#define GLOBAL_FUNCTIONS_GROMACS_HPP

#if 0
namespace gmx_gpu_fmm{

struct fake_aos31_gromacs
{
    typedef real value_type;
    typedef size_t size_type;

    typedef XYZQ<real> Real4;

    const rvec      * xyz;
    const real      * q;
    const size_type n;

    fake_aos31_gromacs(const rvec * xyz, const real * q, const size_type n)
        : xyz(xyz), q(q), n(n)
    { }

    size_type size() const
    {
        return n;
    }

    Real4 operator[] (size_type i) const
    {
        return make_xyzq<real>(xyz[i][0], xyz[i][1], xyz[i][2], q[i]);
    }
};

}//namespace end
#endif

#endif // GLOBAL_FUNCTIONS_GROMACS_HPP
