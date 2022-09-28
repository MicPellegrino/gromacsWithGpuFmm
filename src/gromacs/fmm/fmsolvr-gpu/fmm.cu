#ifndef _BK_fmm_gpu_hpp
#define _BK_fmm_gpu_hpp

#include "fmm.hpp"
#include "cuda_lib.hpp"

namespace gmx_gpu_fmm{

fmm_algorithm::fmm_algorithm():
    p(10), ws(1), depth(3),
    open_boundary_conditions(false),
    dipole_compensation(true),
    printer(this),
    expansion_points(boxes_above_depth(depth + 1)),
    sparse_fmm(false)
    /*
    box_id_map(boxes_above_depth(depth + 1)),
    box_particle_offset(boxes_on_depth(depth) + 1, 0),
    box_particle_offset_host(boxes_on_depth(depth) + 1, 0),
    box_particle_offset_int(boxes_on_depth(depth) + 1, 0),
    particles_per_box_int(boxes_on_depth(depth)),
    q0abc(4, Real4(0.,0.,0.,0.)),
    fake_particles(128, Real4(0.,0.,0.,0.)),
    q0abc_host(4, Real4(0.,0.,0.,0.)),
    fake_particle_size(0),
    fake_particles_host(128, Real4(0.,0.,0.,0.)),
    max_particles_in_box(0),
    initial_box_particle_mem_sizes(boxes_on_depth(depth),0),
    exec_time(), exec_time(),
    lattice_scale(1.0),
    lattice_rescale(1.0)*/ {}

fmm_algorithm::fmm_algorithm(size_t p, ssize_t ws, size_t depth, bool open_boundary, bool compensate_dipole, bool sparse)

    :
      p(p), ws(ws), depth(depth),
      open_boundary_conditions(open_boundary),
      dipole_compensation(compensate_dipole),
      expansion_points(boxes_above_depth(depth + 1)),
      box_id_map(boxes_above_depth(depth + 1)),
      box_particle_offset(boxes_on_depth(depth) + 1, 0),
      box_particle_offset_host(boxes_on_depth(depth) + 1, 0),
      box_particle_offset_int(boxes_on_depth(depth) + 1, 0),
      particles_per_box_int(boxes_on_depth(depth)),
      q0abc(4, Real4(0.,0.,0.,0.)),
      fake_particles(128, Real4(0.,0.,0.,0.)),
      q0abc_host(4, Real4(0.,0.,0.,0.)),
      fake_particle_size(0),
      fake_particles_host(128, Real4(0.,0.,0.,0.)),
      max_particles_in_box(0),
      initial_box_particle_mem_sizes(boxes_on_depth(depth),0),
      printer(this), exec_time(),
      lattice_scale(1.0),
      lattice_rescale(1.0),
      sparse_fmm(sparse)

{
    num_of_all_ops = (3+4*ws)*(3+4*ws)*(3+4*ws);
    num_of_efective_ops = (2+4*ws)*(2+4*ws)*(2+4*ws) - (1+2*ws)*(1+2*ws)*(1+2*ws);
    num_boxes_lowest = boxes_on_depth(depth);
    num_boxes_tree = boxes_above_depth(depth + 1);
    global_offset = boxes_above_depth(depth);
    empty_boxes = 0;
    p1 = p+1;
    p1xx2 = p1*p1;
    p1xp2_2 = (p1*(p1+1))/2;
    pxp1_2 = (p*p1)/2;
    current_stream = 0;
    current_priority_stream = 0;

    ws_dim = (2*ws+1)*2;
    ws_dim2 = ws_dim*ws_dim;
    ws_dim3 = ws_dim2*ws_dim;

    int leastPriority;
    int greatestPriority;

    cudaDeviceGetStreamPriorityRange(&leastPriority, &greatestPriority);

    for (int i=0; i<P2P_STREAMS; ++i)
    {
        CUDA_SAFE_CALL( cudaStreamCreateWithPriority(&streams[i], cudaStreamNonBlocking, greatestPriority));
        CUDA_SAFE_CALL( cudaEventCreateWithFlags(&events[i], cudaEventDisableTiming));
    }

    for (int i=0; i<STREAMS; ++i)
    {
        CUDA_SAFE_CALL( cudaStreamCreateWithPriority(&priority_streams[i], cudaStreamNonBlocking, greatestPriority));
        CUDA_SAFE_CALL( cudaEventCreateWithFlags(&priority_events[i], cudaEventDisableTiming));
    }

    cudaEventCreateWithFlags(&copy_particles_to_host_event, cudaEventBlockingSync);
    cudaEventCreateWithFlags(&ordered_particles_set_event, cudaEventDisableTiming);

    //only exponents in floating point numbers
    //lattice_scale = 1.329227995784916E+36; //2^120
    //lattice_rescale = 7.52316384526264E-37; //1/2^120
    Dipole_compensation = new DipoleCompensation();
}

fmm_algorithm::~fmm_algorithm(){

    //printf("destroy FMM\n");
    cudaDeviceSynchronize();

    for (int i=0; i<P2P_STREAMS; ++i)
    {
        //printf("destroy streams\n");
        CUDA_SAFE_CALL( cudaStreamDestroy( streams[i] ) );
        CUDA_SAFE_CALL(cudaEventDestroy(events[i]));
    }
    for (int i=0; i<STREAMS; ++i)
    {
        CUDA_SAFE_CALL( cudaStreamDestroy( priority_streams[i] ) );
        CUDA_SAFE_CALL(cudaEventDestroy(priority_events[i]));
    }

    cudaEventDestroy(copy_particles_to_host_event);
    cudaEventDestroy(ordered_particles_set_event);

    free_data();

    delete io->result_ptr;
    delete io->result_ptr_host;
    delete Dipole_compensation;
}

void fmm_algorithm::alloc_result(){

    io->result_ptr = new outputadapter_type(io->potential, io->efield);
    io->result_ptr_host = new host_outputadapter_type(io->potential_host, io->efield_host);
}

void fmm_algorithm::update_fmm_boxscale(Real33& abc)
{
    io->box_scale = io->box_size / abc.a.x;
    io->abc       = abc;
    io->abc      *= io->box_scale;

    if(0)
    {
        std::cout<<"x      "<<io->abc.a<<std::endl;
        std::cout<<"y      "<<io->abc.b<<std::endl;
        std::cout<<"z      "<<io->abc.c<<std::endl;
        std::cout<<"corner "<<io->reference_corner<<std::endl;
        std::cout<<"scale  "<<io->box_scale<<std::endl;

        std::cout<<"x      "<<io->half_abc.a<<std::endl;
        std::cout<<"y      "<<io->half_abc.b<<std::endl;
        std::cout<<"z      "<<io->half_abc.c<<std::endl;
    }
}

void fmm_algorithm::init_io(size_t n_particles)
{
    io_unique = std::unique_ptr<input_output<Real> >(new input_output<Real>(n_particles, open_boundary_conditions));
    io = io_unique.get();
    io->reference_corner = Real3(0.,0.,0.);

    size_t size = io->n/32 + num_boxes_lowest;
    offset_map.resize(size);
    block_map.resize(size);
    offset_map_host.resize(size);
    block_map_host.resize(size);

    ordered_particles.resize(io->n);
    ordered_particles_host.resize(io->n);
    orig_ids.resize(io->n);
    fmm_ids.resize(io->n);

    alloc_result();
}

void fmm_algorithm::alloc_exclusions(std::vector<std::vector<int> > &excl){

    if(0)
    {
        for (size_t i = 0; i < excl.size(); ++i)
        {
            printf("i = %d --> ", i);
            for (size_t j = 0; j < excl[i].size(); ++j)
            {
                printf("( j = %d exclusion = %d ) ", j, excl[i][j]);
            }
            printf("\n");
        }
    }

    int excl_n   = 0;
    int mem_size = 0;
    for (size_t i = 0; i < excl.size(); ++i)
    {
        mem_size += excl[i].size();
        excl_n   += excl[i].size() > 0 ? 1 : 0;
    }
    io->excl_n = excl_n;

    if(mem_size == 0)
    {
        alloc_exclusions();
    }
    else
    {
        Device::custom_alloc(exclusions_memory, mem_size*sizeof(int));
        Device::custom_alloc(exclusions,        excl_n*sizeof(int*));
        Device::custom_alloc(exclusions_sizes,  excl_n*sizeof(int));
        int index  = 0;
        for (size_t i = 0; i < excl_n; ++i)
        {
            exclusions_sizes[i] = excl[i].size() - 1;
            exclusions[i] = &exclusions_memory[index];

            for (size_t j = 0; j < excl[i].size(); ++j)
            {
                if(i != excl[i][j])
                {
                    exclusions_memory[index++] = excl[i][j];
                }
            }
        }
    }
    if(0)
    {
        for (size_t i = 0; i < excl.size(); ++i)
        {
            printf("i = %d --> ", i);
            for (size_t j = 0; j < exclusions_sizes[i]; ++j)
            {
                printf("j = %d, exclusion = %d", j, exclusions[i][j]);
            }
            printf("\n");
        }
    }
}

void fmm_algorithm::alloc_exclusions(){

    Device::custom_alloc(exclusions_memory, io->n*sizeof(int));
    Device::custom_alloc(exclusions,        io->n*sizeof(int*));
    Device::custom_alloc(exclusions_sizes,  io->n*sizeof(int));
    Device::custom_alloc(exclusion_pairs,   io->n*sizeof(int));
    for (int i = 0; i < io->n; ++i)
    {
        exclusions[i] = &exclusions_memory[i];
        exclusions_sizes[i] = 0;
        exclusions_memory[i] = i;
    }
}

void fmm_algorithm::distribute_particle_memory(
        const Real33 & abc,
        const Real3 & ref_corner,
        const Real4* input,
        int fmm_depth,
        size_t* mem_size_per_box,
        size_t n)
{
    typedef typename Real3::value_type Real;
    typedef Real4 particle_type;
    typedef uint32_t boxid_type;

    Real33 convert_to_abc = change_of_basis_from_standard_basis(abc);
    Real scale = 1 << fmm_depth;

    boxid_type box_x;
    boxid_type box_y;
    boxid_type box_z;
    boxid_type boxid;
    for (size_t i = 0; i < n; ++i)
    {
        size_t index = i;
        Real3 particle = Real3(input[index]);
        if (particle.x >= abc.a.x)
        {
            particle.x -= abc.a.x;
        }
        if(particle.x < 0.0)
        {
            particle.x += abc.a.x;
        }

        if (particle.y >= abc.b.y)
        {
            particle.y -= abc.b.y;
        }
        if(particle.y < 0.0)
        {
            particle.y += abc.b.y;
        }

        if (particle.z >= abc.c.z)
        {
            particle.z -= abc.c.z;
        }
        if(particle.z < 0.0)
        {
            particle.z += abc.c.z;
        }

        Real3 normalized = convert_to_abc * ( particle - ref_corner );
        
        box_x = normalized.x * scale;
        box_y = normalized.y * scale;
        box_z = normalized.z * scale;
        boxid = make_boxid(box_x, box_y, box_z, fmm_depth);

        mem_size_per_box[boxid]++;
    }
}

void fmm_algorithm::alloc_box_particle_memory(float factor)
{
    Device::devSync();
    distribute_particle_memory(io->abc, io->reference_corner, &io->unordered_particles_host[0], depth, &initial_box_particle_mem_sizes[0], io->n);

    size_t max_mem = 0;
    for(size_t i = 0; i < num_boxes_lowest; ++i)
    {
        if (initial_box_particle_mem_sizes[i] > max_mem)
        {
            max_mem = initial_box_particle_mem_sizes[i];
        }
    }

    if(sparse_fmm)
        max_mem = int(float(max_mem)*factor);
    else
        max_mem = int(float(max_mem + 1)*1.2);

    for(size_t i = 0; i < num_boxes_lowest; ++i)
    {
        initial_box_particle_mem_sizes[i] = max_mem;
    }

    for(size_t i = 0; i < num_boxes_lowest; ++i)
    {
        //printf("box %d, mem = %d mem2 = %d\n",i, ((initial_box_particle_mem_sizes[i])), memsize);
        box[box_id_map[global_offset + i]].alloc_mem((initial_box_particle_mem_sizes[i]));
    }
}

void fmm_algorithm::alloc_data(){

    exec_time.start("ALLOC DATA");
    alloc_data_impl();
    Device::devSync();
    exec_time.stop("ALLOC DATA");
}

void fmm_algorithm::prepare_p2p(int p2p_version){

    exec_time.start("PREPARE P2P");
    prepare_p2p_impl(p2p_version);
    Device::devSync();
    exec_time.stop("PREPARE P2P");
}

void fmm_algorithm::prepare_m2m(){

    exec_time.start("PREPARE M2M");
    prepare_m2m_impl();
    Device::devSync();
    exec_time.stop("PREPARE M2M");
}

void fmm_algorithm::prepare_m2l(){

    exec_time.start("PREPARE M2L");
    prepare_m2l_impl();
    Device::devSync();
    exec_time.stop("PREPARE M2L");
}

void fmm_algorithm::prepare_l2l(){

    exec_time.start("PREPARE L2L");
    prepare_l2l_impl();
    Device::devSync();
    exec_time.stop("PREPARE L2L");
}

void fmm_algorithm::compute_dipole_compensating_particles()
{
    exec_time.start("DIPOLE COMPENSATION");

    if(dipole_compensation)
    {
        Real33 convert_to_abc = change_of_basis_from_standard_basis(io->abc);
        Real33 convert_from_abc = change_of_basis_to_standard_basis(io->abc);

        // compute monopole and (cartesian) dipole for reference point (0,0,0)
        // the dipole is independent of the reference point as long as the monopole is 0
        Real3 dp = get_dipole_from_device(dipole, priority_streams[current_priority_stream]);
        Real3 dp_abc = convert_to_abc * dp;
        Real3 pdp_abc(!open_boundary_conditions ? dp_abc.x : 0.,
                      !open_boundary_conditions ? dp_abc.y : 0.,
                      !open_boundary_conditions ? dp_abc.z : 0.);
        Real3 pdp = convert_from_abc * pdp_abc;
        if (!open_boundary_conditions)
        {
            pdp = dp;
        }

        Dipole_compensation->set_dipole(this, pdp);
        Dipole_compensation->compute_compensation_charges_and_particles(q0abc_host, fake_particles_host, fake_particle_size);
        cudaMemcpyAsync(&q0abc[0], &q0abc_host[0], 4*sizeof(Real4), cudaMemcpyHostToDevice, priority_streams[current_priority_stream]);
        cudaMemcpyAsync(&fake_particles[0], &fake_particles_host[0], fake_particle_size*sizeof(Real4), cudaMemcpyHostToDevice, priority_streams[current_priority_stream]);
    }
    exec_time.stop("DIPOLE COMPENSATION");
}

template <typename Inputdata>
void fmm_algorithm::copy_particles(Inputdata &particles){

    exec_time.start("COPY PARTICLES");
    io->set_particles(particles);
    Device::devSyncDebug();
    exec_time.stop("COPY PARTICLES");
}

void fmm_algorithm::gmx_copy_particles_buffer(REAL* inputparticles, cudaEvent_t* gmx_h2d_ready_event){

    exec_time.start("COPY PARTICLES");
    gmx_copy_particles_buffer_impl(inputparticles, gmx_h2d_ready_event);
    Device::devSyncDebug();
    exec_time.stop("COPY PARTICLES");
}

void fmm_algorithm::distribute_particles(bool gmx_does_neighbor_search){

    exec_time.start("DISTRIBUTE PARTICLES");
    distribute_particles_impl(gmx_does_neighbor_search);
    Device::devSyncDebug();
    exec_time.stop("DISTRIBUTE PARTICLES");
}

void fmm_algorithm::set_expansion_points(){
    for (size_t d = 0; d <= depth; ++d) {
        Real halfstep = reciprocal(Real(size_t(1) << (d + 1)));
        size_t dim = size_t(1) << d;
        size_t offset = boxes_above_depth(d);
        for (size_t i = 0; i < dim; ++i) {
            for (size_t j = 0; j < dim; ++j) {
                for (size_t k = 0; k < dim; ++k) {
                    size_t id = make_boxid(i, j, k, d);
                    expansion_points[offset + id] = io->reference_corner + io->abc * (Real3(2 * i + 1, 2 * j + 1, 2 * k + 1) * halfstep);
                }
            }
        }
    }
}
double fmm_algorithm::energy(Real eps){

    //energy_dump();
    //get rid of this copy here to its own function
    cudaMemcpyAsync(&io->forces_orig_order_host[0], &io->forces_orig_order[0], io->n*sizeof(Real3), cudaMemcpyDeviceToHost, priority_streams[current_priority_stream]);
    Real Ec = energy_impl();
    return Ec*0.5*eps*io->box_scale;

    //return energy_orig_order();
}

void fmm_algorithm::dump_result_in_orig_order(){

    cudaDeviceSynchronize();
    cudaMemcpy(&io->potential_host[0], &io->potential[0], io->n*sizeof(Real), cudaMemcpyDeviceToHost);
    cudaMemcpy(&io->efield_host[0], &io->efield[0], io->n*sizeof(Real3), cudaMemcpyDeviceToHost);

    std::vector<int> fmm_ids(io->n);
    
    for (size_t i = 0; i < io->n; ++i)
    {   
        fmm_ids[orig_ids[i]] = i;   
    }
    Real force_scale = io->box_scale * io->box_scale;
    printf("FORCES>\n");
    for (size_t i = 0; i < io->n; ++i)
    {
        printf("%.8e %.8e %.8e ",io->unordered_particles_host[i].x/io->box_scale,io->unordered_particles_host[i].y/io->box_scale,io->unordered_particles_host[i].z/io->box_scale);
        printf("%.64e %.64e %.64e %.64e\n",io->efield_host[fmm_ids[i]].x*force_scale,io->efield_host[fmm_ids[i]].y*force_scale,io->efield_host[fmm_ids[i]].z*force_scale, io->potential_host[fmm_ids[i]]*io->box_scale);
    }
    printf("FORCES<\n");
}
double fmm_algorithm::energy_orig_order(){

    Real Ec = 0.;
   
    for (size_t i = 0; i < io->n; ++i)
    {
        Real q = io->unordered_particles_host[orig_ids[i]].q;
        Ec += q * io->result_ptr_host->vPhi_ptr[i];
    }
    return Ec*0.5*io->box_scale;
}

void fmm_algorithm::energy_dump(){

    energy_dump_impl(io->result_ptr);
    Device::devSync();
}

void fmm_algorithm::force_dump(int type){

    force_dump_impl(io->result_ptr,type);
    Device::devSync();
}

fmm_algorithm::Real3 fmm_algorithm:: force_l2_norm(){

    Real3 Norm(0.0,0.0,0.0);

    for (size_t i = 0; i < io->n; ++i)
    {
        //Real q = ordered_particles_host[i].q;
        //Real q = 1.0;
        Real3 tmp = io->result_ptr_host->vF[i];
        tmp *= tmp;
        Norm += tmp;
    }

    return Norm.xyz_sqrt()*io->box_scale*io->box_scale;
}

void fmm_algorithm::prepare_data(){

    exec_time.start("PREP DATASTRUCTURES");
    prepare_data_impl();
    Device::devSyncDebug();
    exec_time.stop("PREP DATASTRUCTURES");
}

void fmm_algorithm::p2p_host(int host_device, int version){

    exec_time.start("P2P_HOST");
    p2p_host_impl(host_device, version);
    exec_time.stop("P2P_HOST");
}

void fmm_algorithm::p2p(bool calc_energies, int host_device, int version){

    exec_time.start("P2P");
    if (open_boundary_conditions)
    {
        p2p_impl_open(io->result_ptr);
    }
    else
    {
        p2p_impl_periodic(io->result_ptr, calc_energies, host_device, version);
    }
    Device::devSyncDebug();
    exec_time.stop("P2P");
}

void fmm_algorithm::p2m(){

    exec_time.start("P2M");
    p2m_impl();
    Device::devSyncDebug();
    exec_time.stop("P2M");
}

void fmm_algorithm::m2m(){

    exec_time.start("M2M");
    m2m_impl();
    Device::devSyncDebug();
    exec_time.stop("M2M");
}

void fmm_algorithm::m2l(){

    //std::vector<size_t> chunks = m2l_impl_gpu_cpu_prep();
    exec_time.start("M2L");
    m2l_impl();
    //m2l_impl_gpu_cpu(chunks);
    //m2l_impl_cpu();
    Device::devSyncDebug();
    exec_time.stop("M2L");
}

void fmm_algorithm::lattice(){

    exec_time.start("LATTICE");
    lattice_impl();
    Device::devSyncDebug();
    exec_time.stop("LATTICE");
}

void fmm_algorithm::l2l(){
    exec_time.start("L2L");
    l2l_impl();
    Device::devSyncDebug();
    exec_time.stop("L2L");
}

void fmm_algorithm::forces(Real one4pieps0){
    exec_time.start("FORCES");
    forces_impl(io->result_ptr, one4pieps0);
    Device::devSyncDebug();
    exec_time.stop("FORCES");
}

void fmm_algorithm::get_forces_in_orig_order(float3* gmx_forces, cudaStream_t* gmx_sync_stream){

    get_forces_in_orig_order_impl(gmx_forces, gmx_sync_stream);
}

std::pair<double, double> fmm_algorithm::get_error(fmm_algorithm::CoeffMatrix **sol, fmm_algorithm::CoeffMatrix **ref){

    double relative_real_error = 0.0;
    double relative_imag_error = 0.0;

    double real_norm = 0;
    double imag_norm = 0;
    double real_error = 0;
    double imag_error = 0;

    double real = 0;
    double imag = 0;

    for (size_t d = 0; d <= depth; ++d) {

        size_t depth_offset = boxes_above_depth(d);
        ssize_t dim = ssize_t(1) << d;
        for (ssize_t i = 0; i < dim; ++i) {
            for (ssize_t j = 0; j < dim; ++j) {
                for (ssize_t k = 0; k < dim; ++k) {
                    size_t id = depth_offset + make_boxid(i, j, k, d);

                    for (int l = 0; l <= p; ++l) {
                        for (int m = 0; m <= l; ++m) {

                            real =((*ref[id])(l,m)).real()-((*sol[id])(l,m)).real();
                            imag =((*ref[id])(l,m)).imag()-((*sol[id])(l,m)).imag();

                            real_norm +=((*ref[id])(l,m)).real()*((*ref[id])(l,m)).real();
                            imag_norm +=((*ref[id])(l,m)).imag()*((*ref[id])(l,m)).imag();

                            real_error +=real*real;
                            imag_error +=imag*imag;

                        }
                    }
                }
            }
        }
    }

    if(real_norm != 0.0)
        relative_real_error = std::pow(real_error/real_norm,0.5);
    if(imag_norm != 0.0)
        relative_imag_error = std::pow(imag_error/imag_norm,0.5);

    std::pair<double,double> errors(relative_real_error,relative_imag_error);

    return errors;
}

template void fmm_algorithm::copy_particles<fake_AoS<fmm_algorithm::Real>>(fake_AoS<fmm_algorithm::Real> & particle_data);

template void fmm_algorithm::copy_particles<fake_AoS_GMX<fmm_algorithm::Real>>(fake_AoS_GMX<fmm_algorithm::Real> & particle_data);
template void fmm_algorithm::copy_particles<fake_AoS_GMX<double, float>>(fake_AoS_GMX<double, float> & particle_data);
}//namespace end

#endif
