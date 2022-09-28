#ifndef _BK_FMM_dipole_compensation_hpp_
#define _BK_FMM_dipole_compensation_hpp_

#include <algorithm>

#include "linear_algebra.hpp"
#include "particle2multipole.hpp"
#include "particle2local.hpp"
#include "cuda_DIPOLE_COMPENSATION.hpp"

namespace gmx_gpu_fmm{

namespace fmsolvr {

template <typename SimulationBox, typename OutputContainer>
class fake_particle_generator
{
    typedef typename OutputContainer::value_type particle_type;
    typedef typename particle_type::value_type Real;
    typedef typename SimulationBox::Real3 Real3;

    const SimulationBox & simulationbox_;
    OutputContainer & output_;
    Real3 &root_;

public:
    size_t index_;

    fake_particle_generator(const SimulationBox & simbox, OutputContainer & output, Real3 & root) :
        simulationbox_(simbox), output_(output), root_(root), index_(0)
    { }

    void generate(size_t ws, Real q0, Real qa, Real qb, Real qc)
    {
        gen(q0, -ws, -ws, -ws, -ws, -ws, -ws);
        gen(qa, ws+1, ws+1, -ws, +ws, -ws, +ws);
        gen(qb, -ws, +ws, ws+1, ws+1, -ws, +ws);
        gen(qc, -ws, +ws, -ws, +ws, ws+1, ws+1);
        gen(q0 + qa, -ws+1, +ws, -ws, -ws, -ws, -ws);
        gen(q0 + qb, -ws, -ws, -ws+1, +ws, -ws, -ws);
        gen(q0 + qc, -ws, -ws, -ws, -ws, -ws+1, +ws);
        gen(q0 + qa + qb, -ws+1, +ws, -ws+1, +ws, -ws, -ws);
        gen(q0 + qa + qc, -ws+1, +ws, -ws, -ws, -ws+1, +ws);
        gen(q0 + qb + qc, -ws, -ws, -ws+1, +ws, -ws+1, +ws);
    }

private:
    void gen(Real q, int alo, int ahi, int blo, int bhi, int clo, int chi)
    {
        Real3 shift = simulationbox_.reference_corner;

        for (int i = alo; i <= ahi; ++i)
            for (int j = blo; j <= bhi; ++j)
                for (int k = clo; k <= chi; ++k) {
                    Real3 ijk(i, j, k);
                    Real3 pos = simulationbox_.abc * ijk + shift;
                    output_[index_++] = make_xyzq(Real3(pos.x, pos.y, pos.z) - root_ , q);
                }
    }
};

template <typename SimulationBox, typename OutputContainer, typename Real, typename Real3>
void generate_dipole_compensating_particles(OutputContainer & output, const SimulationBox & simbox, size_t ws, Real q0, Real qa, Real qb, Real qc, size_t &size, Real3 root)
{
    fake_particle_generator<SimulationBox, OutputContainer> gfp(simbox, output, root);
    gfp.generate(ws, q0, qa, qb, qc);
    size = gfp.index_;
}

template <typename Real, typename Real3, typename Real33>
void compute_dipole_compensating_charges(const Real3 &pdp, const Real33 &abc, Real &q0, Real &qa, Real & qb, Real & qc)
{
    linear_algebra::SqMatrix<Real> A(4);
    Real3 origin(0., 0., 0.);
    A[0][0] = A[0][1] = A[0][2] = A[0][3] = 1;
    A[1][0] = origin.x;
    A[1][1] = abc.a.x;
    A[1][2] = abc.b.x;
    A[1][3] = abc.c.x;
    A[2][0] = origin.y;
    A[2][1] = abc.a.y;
    A[2][2] = abc.b.y;
    A[2][3] = abc.c.y;
    A[3][0] = origin.z;
    A[3][1] = abc.a.z;
    A[3][2] = abc.b.z;
    A[3][3] = abc.c.z;
    Real dipole[4] = { 0, -pdp.x, -pdp.y, -pdp.z };
    Real fiqtive[4];
    if (!linear_algebra::solve(4, A, fiqtive, dipole))
            throw "eqn system does not have a unique solution";
    q0 = fiqtive[0];
    qa = fiqtive[1];
    qb = fiqtive[2];
    qc = fiqtive[3];
    //linear_algebra::dump(4, A, fiqtive, dipole);
}

template <typename FMMHandle>
class dipole_compensation
{
    typedef typename FMMHandle::Real33 Real33;
    typedef typename FMMHandle::Real4 Real4;
    typedef typename FMMHandle::Real3 Real3;
    typedef typename FMMHandle::Real Real;
    typedef size_t  size_type;
    typedef typename FMMHandle::io_type simulationbox_type;

public:
    dipole_compensation()
    {}

    void set_dipole(FMMHandle* fmm_handle, Real3& dipole)
    {
        fmm_handle_ = fmm_handle;
        dipole_ = dipole;
    }

    template <typename particle_vector>
    void compute_compensation_charges_and_particles(particle_vector& q0abc, particle_vector& fake_particles, size_t &fake_particles_size)
    {

        compute_dipole_compensating_charges(dipole_, fmm_handle_->io->abc, q0_, qa_, qb_, qc_);
        Real33 abc = fmm_handle_->io->abc;
        Real3 origin(0., 0., 0.);
        Real3 m = (abc.a + abc.b + abc.c) * 0.5;

        q0abc[0] = make_xyzq(origin - m, q0_);
        q0abc[1] = make_xyzq(abc.a - m , qa_);
        q0abc[2] = make_xyzq(abc.b - m, qb_);
        q0abc[3] = make_xyzq(abc.c - m, qc_);
        Real3 expn_root = m;

        //std::cout<<q0abc[0]<<q0abc[1]<<q0abc[2]<<q0abc[3]<<std::endl;

        generate_dipole_compensating_particles(fake_particles, *fmm_handle_->io, fmm_handle_->ws, q0_, qa_, qb_, qc_, fake_particles_size, expn_root);

    }

    Real q0() const
    { return q0_; }

    Real qa() const
    { return qa_; }

    Real qb() const
    { return qb_; }

    Real qc() const
    { return qc_; }

private:
    FMMHandle* fmm_handle_;
    Real3 dipole_;
    Real q0_, qa_, qb_, qc_;
};

}  // namespace fmsolvr

}//namespace end

#endif
// vim: et:ts=4:sw=4

