#include "fmm.hpp"

namespace gmx_gpu_fmm{

void fmm_algorithm::alloc_data_impl(){

    //P2P memory
    Device::custom_alloc(box_particle_source_starts,  (num_boxes_lowest*27)*sizeof(size_t));
    Device::custom_alloc(box_particle_source_ends, (num_boxes_lowest*27)*sizeof(size_t));
    Device::custom_alloc(box_particle_source_shifts, (num_boxes_lowest*27)*sizeof(Real3));
    Device::custom_alloc(box_particle_source_offsets,  (num_boxes_lowest*27)*sizeof(size_t));

    //A operator memory
    Device::custom_alloc(box_a_operators,  (num_boxes_tree-1)*sizeof(M2M_Operator*));
    Device::custom_alloc(box_a_targets, (num_boxes_tree-1)*sizeof(CoeffMatrix*));
    Device::custom_alloc(box_a_target_ids, (num_boxes_tree-1)*sizeof(size_t));
#ifdef M2M_SOA_OPTIMIZATION
    Device::custom_alloc(box_a_targets_SoA, p1xp2_2*sizeof(value_type**));

    for (size_t i = 0; i<p1xp2_2;++i){
        Device::custom_alloc(box_a_targets_SoA[i], (num_boxes_tree-1)*sizeof(value_type*));
    }
#endif
    //B operator memory
    Device::custom_alloc(box_b_operators, (num_boxes_tree-1)*num_of_efective_ops*sizeof(M2L_Operator*));
    Device::custom_alloc(box_b_targets, (num_boxes_tree-1)*num_of_efective_ops*sizeof(CoeffMatrix*));
    Device::custom_alloc(box_b_target_ids, (num_boxes_tree-1)*num_of_efective_ops*sizeof(size_t));
#ifdef M2L_SOA_OPTIMIZATION
    Device::custom_alloc(box_b_targets_SoA, p1xp2_2*sizeof(value_type**));

    for (size_t i = 0; i<p1xp2_2;++i){
        Device::custom_alloc(box_b_targets_SoA[i], (num_boxes_tree-1)*num_of_efective_ops*sizeof(value_type*));
    }
#endif
    //C operator memory
    Device::custom_alloc(box_c_operators, (num_boxes_tree-1)*8*sizeof(L2L_Operator*));
    Device::custom_alloc(box_c_targets, (num_boxes_tree-1)*8*sizeof(CoeffMatrix*));
    Device::custom_alloc(box_c_target_ids, (num_boxes_tree-1)*8*sizeof(size_t));
#ifdef L2L_SOA_OPTIMIZATION
    Device::custom_alloc(box_c_targets_SoA, p1xp2_2*sizeof(value_type**));

    for (size_t i = 0; i<p1xp2_2;++i){
        Device::custom_alloc(box_c_targets_SoA[i], (num_boxes_tree-1)*8*sizeof(value_type*));
    }
#endif
    Device::custom_alloc(box_operators_tmp, num_of_efective_ops*sizeof(M2L_Operator*));
    Device::custom_alloc(box_targets_tmp, num_of_efective_ops*sizeof(CoeffMatrix*));
    Device::custom_alloc(box_target_ids_tmp, num_of_efective_ops*sizeof(size_t));

    Device::custom_alloc(m2m_operators, depth*sizeof(M2M_Operator*));
    Device::custom_alloc(m2l_operators, depth*sizeof(M2L_Operator*));
    Device::custom_alloc(l2l_operators, depth*sizeof(L2L_Operator*));

    Device::custom_alloc(dummy_m2l_operators, depth*sizeof(M2L_Operator*));

    for (size_t d_c = depth; d_c > 0; --d_c)
    {
        size_t d_p = d_c - 1;
        size_t d_delta = d_c - d_p;
        size_t num_operators = boxes_on_depth(d_delta);
        Device::custom_alloc(m2m_operators[d_p], num_operators*sizeof(M2M_Operator*));
        Device::custom_alloc(l2l_operators[d_p], num_operators*sizeof(L2L_Operator*));
        Device::custom_alloc(m2l_operators[d_p], num_of_all_ops*sizeof(M2L_Operator*));
    }

    Device::custom_alloc(omega, num_boxes_tree*sizeof(CoeffMatrix*));
    Device::custom_alloc(mu, num_boxes_tree*sizeof(CoeffMatrix*));

    omega_memory = new CoeffMatrixMemory(p,num_boxes_tree);
    omegav       = new CoeffMatrix[num_boxes_tree];
    mu_memory    = new CoeffMatrixMemory(p,num_boxes_tree);
    muv          = new CoeffMatrix[num_boxes_tree];

    for (size_t i=0;i<num_boxes_tree;i++)
    {
        //omega[i] = new CoeffMatrix(p);
        //mu[i] = new CoeffMatrix(p);
        omegav[i].reinit(p,omega_memory->get_raw_pointer(i));
        muv[i].reinit(p,mu_memory->get_raw_pointer(i));

        omega[i] = &omegav[i];
        mu[i] = &muv[i];
    }

    //host
    Host::custom_alloc(omega_host, num_boxes_tree*sizeof(CoeffMatrix_Host*));
    Host::custom_alloc(mu_host, num_boxes_tree*sizeof(CoeffMatrix_Host*));
    //omega_memory_host = new CoeffMatrixMemory_Host(p,num_boxes_tree);
    //omegav_host       = new CoeffMatrix_Host[num_boxes_tree];
    //mu_memory_host    = new CoeffMatrixMemory_Host(p,num_boxes_tree);
    //muv_host          = new CoeffMatrix_Host[num_boxes_tree];

    for (size_t i=0;i<num_boxes_tree;i++)
    {
        omega_host[i] = new CoeffMatrix_Host(p);
        mu_host[i] = new CoeffMatrix_Host(p);
        //omegav_host[i].reinit(p,omega_memory_host->get_raw_pointer(i));
        //muv_host[i]->reinit(p,mu_memory_host->get_raw_pointer(i));

        //omega_host[i] = &omegav_host[i];
        //mu_host[i] = &muv_host[i];
    }


    dipole_host = Real3(0.,0.,0.);
    DeviceOnly::custom_alloc(dipole, sizeof(Real3));

    Lattice = new CoeffMatrix(2 * p);
    LatticeD = new CoeffMatrixD(2 * p);
    if (!open_boundary_conditions)
    {
        typedef typename CoeffMatrixD::complex_type complex_type;
        typedef ABC<XYZ<double> > Real33D;
        XYZ<double> e1((double)io->abc.a.x,(double)io->abc.a.y,(double)io->abc.a.z);
        XYZ<double> e2((double)io->abc.b.x,(double)io->abc.b.y,(double)io->abc.b.z);
        XYZ<double> e3((double)io->abc.c.x,(double)io->abc.c.y,(double)io->abc.c.z);
        Real33D abcD(e1, e2, e3);
        size_t its = LatticeOperator3D<Real33D, CoeffMatrixD>(abcD, ws, *LatticeD, p, 42);
        Real imag_scale = 0.0;
        if(lattice_scale != 1.0)
            imag_scale = 1.0;
        *LatticeD *= complex_type(lattice_scale, imag_scale);
        //printf("LATTICE \n");
        //dump(*LatticeD,2*p);
        //printf("Lattice max %e\n", LatticeD->max());
        //printf("Lattice min %e\n", LatticeD->min());
        //Lattice->populate_lower();
        Lattice->recast(*LatticeD);
    }

    box = new Box[num_boxes_tree + 1];

    if(depth == 0)
    {
        box[0].init(0,0);
        box[0].set_offset_mem(0, box_particle_source_offsets, box_particle_source_shifts);
    }
    else
    {
        box[0].init(0,0);
    }
    box[0].set_omega(omega[0]);
    box[0].set_mu(mu[0]);
    box[0].set_c_mem(0 * 8, box_c_operators, box_c_targets, box_c_target_ids);
    box_id_map[0] = 0;
    size_t box_id = 1;
    size_t deepest_level_index = 0;

    //should be computed
    //size_t particles_per_box = io->n;

    for (size_t d = 1; d <= depth; ++d)
    {
        ssize_t dim = ssize_t(1) << d;
        size_t depth_offset = boxes_above_depth(d);
        for (ssize_t z = 0; z < 2; ++z)
        {
            for (ssize_t y = 0; y < 2; ++y)
            {
                for (ssize_t x = 0; x < 2; ++x)
                {
                    for (ssize_t i = z; i < dim; i+=2)
                    {
                        for (ssize_t j = y; j < dim; j+=2)
                        {
                            for (ssize_t k = x; k < dim; k+=2)
                            {
                                size_t omega_id = depth_offset + make_boxid(i, j, k, d);

                                if(d<depth)
                                {
                                    box[box_id].init(omega_id, d);
                                }
                                else
                                {
                                    box[box_id].init(omega_id, d);
                                    size_t particles_memory_index = deepest_level_index*27;
                                    box[box_id].set_offset_mem(particles_memory_index, box_particle_source_offsets, box_particle_source_shifts);
                                    deepest_level_index++;
                                }
                                box[box_id].set_a_mem((box_id-1), box_a_operators, box_a_targets, box_a_target_ids);
                                box[box_id].set_b_mem((box_id-1) * num_of_efective_ops, box_b_operators, box_b_targets, box_b_target_ids);
                                box[box_id].set_c_mem(box_id * 8, box_c_operators, box_c_targets, box_c_target_ids);
                                box[box_id].set_omega(omega[omega_id]);
                                box[box_id].set_mu(mu[omega_id]);
                                box_id_map[omega_id] = box_id;
                                box_id++;
                            }
                        }
                    }
                }
            }
        }
    }

    typedef typename CoeffMatrix::value_type value_type;
    typedef typename Device::allocator::template rebind<value_type>::other alloc_complex;
    omegaSoA = new CoeffMatrixSoA(p);
    muSoA = new CoeffMatrixSoA(p);

    for(ssize_t l = 0; l <= (ssize_t)p; ++l)
    {
       for(ssize_t m = 0; m <= l; ++m)
       {
           Device::custom_alloc(omegaSoA->getSoA(l,m), num_boxes_tree*sizeof(value_type));
           Device::custom_alloc(muSoA->getSoA(l,m), num_boxes_tree*sizeof(value_type));
           for(int i = 0; i < num_boxes_tree; i++)
           {
               omegaSoA->getSoA(l,m)[i] = value_type(0.,0.);
               muSoA->getSoA(l,m)[i] = value_type(0.,0.);
           }
       }
    }

    //needed for forces
    Device::custom_alloc(dmux, num_boxes_lowest*sizeof(CoeffMatrix*));
    Device::custom_alloc(dmuy, num_boxes_lowest*sizeof(CoeffMatrix*));
    Device::custom_alloc(dmuz, num_boxes_lowest*sizeof(CoeffMatrix*));

    dmux_memory = new CoeffMatrixMemory(p, num_boxes_lowest);
    dmuy_memory = new CoeffMatrixMemory(p, num_boxes_lowest);
    dmuz_memory = new CoeffMatrixMemory(p, num_boxes_lowest);

    dmuxv       = new CoeffMatrix[num_boxes_lowest];
    dmuyv       = new CoeffMatrix[num_boxes_lowest];
    dmuzv       = new CoeffMatrix[num_boxes_lowest];

    for (size_t i=0;i<num_boxes_lowest;i++)
    {
        //dmux[i] = new CoeffMatrix(p);
        //dmuy[i] = new CoeffMatrix(p);
        //dmuz[i] = new CoeffMatrix(p);

        dmuxv[i].reinit(p, dmux_memory->get_raw_pointer(i));
        dmuyv[i].reinit(p, dmuy_memory->get_raw_pointer(i));
        dmuzv[i].reinit(p, dmuz_memory->get_raw_pointer(i));

        dmux[i] = &dmuxv[i];
        dmuy[i] = &dmuyv[i];
        dmuz[i] = &dmuzv[i];
    }

    //Device::custom_alloc(sparse_interactions_list, num_boxes_tree*189*sizeof(int2));
    //Device::custom_alloc(sparse_interactions_list_size,1*sizeof(int));

    //memory for energy chunks reduced on the device
    //device starts 512 blocks of 32 stripes = 16384
    //Ec will work for maximal 268435456 particles
    DeviceOnly::custom_alloc(Ec,65536*sizeof(Real));
}

void fmm_algorithm::free_data(){

    Device::custom_free(exclusions_memory);
    Device::custom_free(exclusions);
    Device::custom_free(exclusions_sizes);

    //P2P memory
    Device::custom_free(box_particle_source_starts);
    Device::custom_free(box_particle_source_ends);
    Device::custom_free(box_particle_source_shifts);
    Device::custom_free(box_particle_source_offsets);

    if(depth > 0)
    {
#ifdef M2M_SOA_OPTIMIZATION
        for (size_t i = 0; i<p1xp2_2;++i){
            Device::custom_free(box_a_targets_SoA[i]);
        }
#endif
        Device::custom_free(box_a_operators);
        Device::custom_free(box_a_targets);
        Device::custom_free(box_a_target_ids);

#ifdef M2L_SOA_OPTIMIZATION
        //B operator memory
        for (size_t i = 0; i<p1xp2_2;++i){
            Device::custom_free(box_b_targets_SoA[i]);
        }
#endif

        Device::custom_free(box_b_operators);
        Device::custom_free(box_b_targets);
        Device::custom_free(box_b_target_ids);
#ifdef L2L_SOA_OPTIMIZATION
        //C operator memory
        for (size_t i = 0; i<p1xp2_2;++i){
            Device::custom_free(box_c_targets_SoA[i]);
        }
#endif
        Device::custom_free(box_c_operators);
        Device::custom_free(box_c_targets);
        Device::custom_free(box_c_target_ids);
    }

#ifdef M2M_SOA_OPTIMIZATION
    Device::custom_free(box_a_targets_SoA);
#endif

#ifdef M2L_SOA_OPTIMIZATION
    Device::custom_free(box_b_targets_SoA);
#endif

#ifdef L2L_SOA_OPTIMIZATION
    Device::custom_free(box_c_targets_SoA);
#endif

    Device::custom_free(box_operators_tmp);
    Device::custom_free(box_targets_tmp);
    Device::custom_free(box_target_ids_tmp);

    delete omega_memory;
    delete mu_memory;

    for (size_t i=0;i<num_boxes_tree;i++)
    {
        //delete omega[i];
        //delete mu[i];
    }

    Device::custom_free(omega);
    Device::custom_free(mu);

    delete Lattice;
    delete LatticeD;

    for(size_t i = 0; i < num_boxes_lowest; ++i)
    {
        box[global_offset + i].free_mem();
    }

    delete[] box;

    if(depth > 0)
    {
        for (size_t d_c = depth; d_c > 0; --d_c)
        {
            size_t d_p = d_c - 1;
            for (int idx_op = 0; idx_op < 8; ++idx_op)
            {
                delete m2m_operators[d_p][idx_op];
                delete l2l_operators[d_p][idx_op];
            }

            for (int idx_op = 0; idx_op < num_of_all_ops; ++idx_op)
            {
                //printf("depth %d index %d pointer %p\n", d_p, idx_op, m2l_operators[d_p][idx_op]);
                delete m2l_operators[d_p][idx_op];
            }
        }
        Device::custom_free(dummy_m2l_operators);

        for (size_t d_c = depth; d_c > 0; --d_c)
        {
            size_t d_p = d_c - 1;
            Device::custom_free(m2m_operators[d_p]);
            Device::custom_free(l2l_operators[d_p]);
            Device::custom_free(m2l_operators[d_p]);
        }

        Device::custom_free(m2m_operators);
        Device::custom_free(m2l_operators);
        Device::custom_free(l2l_operators);
    }

    delete non_periodic_mu_dummy;
    delete non_periodic_m2l_dummy;

    delete[] bitset;
    Device::custom_free(bitset_mem);

    //host
    for (size_t i=0;i<num_boxes_tree;i++)
    {
        delete omega_host[i];
        delete mu_host[i];
    }

    //delete omega_memory_host;
    //delete mu_memory_host;

    Host::custom_free(omega_host);
    Host::custom_free(mu_host);

    for(ssize_t l = 0; l <= (ssize_t)p; ++l)
    {
       for(ssize_t m = 0; m <= l; ++m)
       {
           Device::custom_free(omegaSoA->getSoA(l,m));
           Device::custom_free(muSoA->getSoA(l,m));
       }
    }

    delete omegaSoA;
    delete muSoA;

    //delete[] dmuxv;
    //delete[] dmuyv;
    //delete[] dmuzv;

    for (size_t i=0;i<num_boxes_lowest;i++)
    {        
        //delete dmux[i];
        //delete dmuy[i];
        //delete dmuz[i];
    }

    delete dmux_memory;
    delete dmuy_memory;
    delete dmuz_memory;

    Device::custom_free(dmux);
    Device::custom_free(dmuy);
    Device::custom_free(dmuz);

if(p2p_version == 1)
{
    for(int i = 0; i < 13; ++i)
    {
        Device::custom_free(p2p_particles_box_pairs[i]);
        Device::custom_free(p2p_particles_periodic_shifts[i]);

        Host::custom_free(p2p_particles_box_pairs_host[i]);
        Host::custom_free(p2p_particles_periodic_shifts_host[i]);
    }

    Device::custom_free(p2p_particles_box_pairs);
    Device::custom_free(p2p_particles_periodic_shifts);

    Host::custom_free(p2p_particles_box_pairs_host);
    Host::custom_free(p2p_particles_periodic_shifts_host);

    Host::custom_free(p2p_particles_host);
    Host::custom_free(p2p_results_host);

    for(size_t i = 0; i < num_boxes_lowest; ++i)
    {
        Device::custom_free(p2p_particles[i]);
        Device::custom_free(p2p_results[i]);
    }

    Device::custom_free(p2p_particles);
    Device::custom_free(p2p_results);
    Device::custom_free(p2p_particles_sizes_div8);
    Device::custom_free(p2p_rounded_particles_sizes);
}
    Device::custom_free(Ec);
}

}//namespace end


