#include "input_output.hpp"
#include "ioadapter.hpp"
#include <omp.h>

namespace gmx_gpu_fmm{

template <typename T>
input_output<T>::input_output(const size_t n_in, bool open_boundary_conditions): n(n_in),

    periodic(!open_boundary_conditions),

    unordered_particles(n_in, Real4()),
    potential(n_in, 0.),
    efield(n_in, Real3(0., 0., 0.)),

    unordered_particles_host(n_in, Real4()),
    potential_host(n_in, 0.),
    forces_orig_order(n_in,Real3()),
    forces_orig_order_host(n_in,Real3()),
    potential_orig_order(n_in, 0.0),
    potential_orig_order_host(n_in, 0.0),
    efield_host(n_in, Real3(0., 0., 0.)),
    potential_host_p2p(n_in, 0.),
    efield_host_p2p(n_in, Real3(0., 0., 0.)),
    dipole(0.,0.,0.),
    box_scale(1.0),
    box_size(32.0)
{
    //printf("box size %e\n", box_size );
    abc.a = Real3(box_size, 0.0, 0.0);
    abc.b = Real3(0.0, box_size, 0.0);
    abc.c = Real3(0.0, 0.0, box_size);

    half_abc = abc.half();
}

template <typename T>
template<typename inputdata>
void input_output<T>::input_output::set_particles(inputdata &particle_data)
{
    //omp_set_num_threads(6);
#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {

        unordered_particles_host[i].x = particle_data.x(i) * box_scale;
        unordered_particles_host[i].y = particle_data.y(i) * box_scale;
        unordered_particles_host[i].z = particle_data.z(i) * box_scale;
        unordered_particles_host[i].q = particle_data.q(i);
        //printf("%e\n", particle_data.q(i));
    }

    if(0)
    {
        for (int i = 0; i < n; ++i)
        {
            std::cout<<particle_data.x(i)<<" "<<particle_data.y(i)<<" "<<particle_data.z(i)<<" --- ";
            std::cout<<unordered_particles_host[i].x<<" "<<unordered_particles_host[i].y<<" "<<unordered_particles_host[i].z<<std::endl;
        }

    }
}

template <typename T>
bool input_output<T>::periodic_a()
{
    return periodic;
}

template <typename T>
bool input_output<T>::periodic_b()
{
    return periodic;
}

template <typename T>
bool input_output<T>::periodic_c()
{
    return periodic;
}

template class input_output<float>;
template class input_output<double>;

typedef fake_AoS<float> fake_aosF;
typedef fake_AoS<double> fake_aosD;

template void input_output<float>::set_particles<fake_aosF>(fake_aosF &particle_dat);
template void input_output<double>::set_particles<fake_aosD>(fake_aosD &particle_data);

typedef fake_AoS_GMX<float> fake_aos_GMX_F;
typedef fake_AoS_GMX<double> fake_aos_GMX_D;
typedef fake_AoS_GMX<double, float> fake_aos_GMX_DF;

template void input_output<float>::set_particles<fake_aos_GMX_F>(fake_aos_GMX_F &particle_data);
template void input_output<float>::set_particles<fake_aos_GMX_DF>(fake_aos_GMX_DF &particle_data);
template void input_output<double>::set_particles<fake_aos_GMX_D>(fake_aos_GMX_D &particle_data);
template void input_output<double>::set_particles<fake_aos_GMX_DF>(fake_aos_GMX_DF &particle_data);

}//namespace end
