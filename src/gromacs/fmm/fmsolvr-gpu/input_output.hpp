#ifndef INPUT_OUTPUT_HPP
#define INPUT_OUTPUT_HPP

#include <cassert>
#include <vector>
#include "data_type.hpp"
#include "cuda_alloc.hpp"
#include "xyz.hpp"
#include "xyzq.hpp"
#include "abc.hpp"
#include "ioadapter.hpp"
#include "global_functions.hpp"
#include "cudavector.hpp"
#include "architecture.hpp"

namespace gmx_gpu_fmm{

template <typename T>
class input_output{

public:

    typedef T Real;

    typedef XYZ<Real>  Real3;
    typedef XYZQ<Real> Real4;
    typedef ABC<Real3> Real33;

    typedef cudavector<Real4, cuda_device_allocator<Real4> >                   particle_vector_type;
    typedef cudavector<Real3, cuda_device_allocator<Real3> >                   field_vector_type;
    //typedef cudavector<Real3, cuda_allocator<Real3> >                          force_vector_type;
    typedef cudavector<Real,  cuda_device_allocator<Real> >                    potential_vector_type;

    typedef std::vector<Real4, cuda_host_allocator<Real4> >                    host_particle_vector_type;   
    typedef std::vector<Real3, cuda_host_allocator<Real3> >                    host_field_vector_type;
    typedef std::vector<Real, cuda_host_allocator<Real> >                      host_potential_vector_type;

    //typedef inputadapter<particle_vector_type>                                 inputadapter_type;
    //typedef inputadapter<host_particle_vector_type>                            host_inputadapter_type;
    typedef outputadapter<potential_vector_type, field_vector_type, Device<Real> >         outputadapter_type;
    typedef outputadapter<host_potential_vector_type, host_field_vector_type, Host<Real>>  host_outputadapter_type;

    potential_vector_type potential;
    field_vector_type efield;
    outputadapter_type* result_ptr;

    host_potential_vector_type potential_host;
    host_field_vector_type efield_host;
    host_outputadapter_type* result_ptr_host;

    host_potential_vector_type potential_host_p2p;
    host_field_vector_type efield_host_p2p;

    particle_vector_type unordered_particles;
    host_particle_vector_type unordered_particles_host;

    field_vector_type forces_orig_order;
    host_field_vector_type forces_orig_order_host;

    potential_vector_type potential_orig_order;
    host_potential_vector_type potential_orig_order_host;

    size_t excl_n;
    size_t n;
    Real33 abc;
    Real33 half_abc;
    Real3 reference_corner;
    bool periodic;
    Real3 dipole;

    Real box_scale;
    Real box_size;

    input_output(const size_t n_in, bool open_boundary_conditions);

    template<typename inputdata>
    void set_particles(inputdata &particle_data);

    bool periodic_a();

    bool periodic_b();

    bool periodic_c();
};

}//namespace end

#endif // INPUT_OUTPUT_HPP
