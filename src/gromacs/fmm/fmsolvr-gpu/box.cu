#include "box.hpp"

namespace gmx_gpu_fmm
{

template <typename Real, typename arch>
Box<Real, arch>::Box()
{
    ptcl_index   = 0;
    active       = 0;
}

template <typename Real, typename arch>
DEVICE
void Box<Real, arch>::set_orig_index(size_t orig_index)
{
    orig_ptcl_ids[atomicAdd(&ptcl_index, 1)] = orig_index;
    active = 1;
}

template <typename Real, typename arch>
void Box<Real, arch>::alloc_mem(size_t n)
{
    arch::custom_alloc(orig_ptcl_ids, n * sizeof(size_t));
    //cudaMemPrefetchAsync(orig_ptcl_ids, n * sizeof(size_t), 0, 0);
}

template <typename Real, typename arch>
void Box<Real, arch>::free_mem()
{
    arch::custom_free(orig_ptcl_ids);
}

template <typename Real, typename arch>
void Box<Real, arch>::init(size_t id_, size_t d)
{
    depth        = d;
    id           = id_;
    depth_offset = boxes_above_depth(d);
}

template <typename Real, typename arch>
void Box<Real, arch>::set_offset_mem(size_t mem_index, size_t* offset_mem, Real3* shifts_mem)
{
    particle_offset_ids      = &offset_mem[mem_index];
    particle_periodic_shifts = &shifts_mem[mem_index];
}

template <typename Real, typename arch>
void Box<Real, arch>::set_a_mem(size_t op_mem_index, M2M_Operator** o_mem, CoeffMatrix** t_mem, size_t* id_mem)
{
    a_operators  = &o_mem[op_mem_index];
    a_targets    = &t_mem[op_mem_index];
    a_target_ids = &id_mem[op_mem_index];
}

template <typename Real, typename arch>
void Box<Real, arch>::set_b_mem(size_t op_mem_index, M2L_Operator** o_mem, CoeffMatrix** t_mem, size_t* id_mem)
{
    b_operators = &o_mem[op_mem_index];
    b_targets   = &t_mem[op_mem_index];
    b_target_ids = &id_mem[op_mem_index];
}

template <typename Real, typename arch>
void Box<Real, arch>::set_c_mem(size_t op_mem_index, L2L_Operator** o_mem, CoeffMatrix** t_mem, size_t *id_mem)
{
    c_operators = &o_mem[op_mem_index];
    c_targets   = &t_mem[op_mem_index];
    c_target_ids = &id_mem[op_mem_index];
}

template <typename Real, typename arch>
void Box<Real, arch>::set_permutations(std::vector<std::vector<size_t> > &permutations)
{
    op_perm = permutations;
}

template <typename Real, typename arch>
void Box<Real, arch>::permute_ops(M2L_Operator** b_operators_temp, CoeffMatrix **targets_temp, size_t* target_ids_temp, size_t num_of_efective_ops)
{
    b_operators_tmp = b_operators_temp;
    targets_tmp     = targets_temp;
    target_ids_tmp  = target_ids_temp;

    for (size_t i = 0; i < num_of_efective_ops; ++i)
    {
        b_operators_tmp[i] = b_operators[i];
        targets_tmp[i]     = b_targets[i];
        target_ids_tmp[i]  = b_target_ids[i];
    }

    size_t index = 0;
    for (size_t type = 0; type < 4; ++type)
    {
        size_t number_in_groups = 0;
        for (size_t i = 0; i < op_perm[type].size(); ++i)
        {
            if (type == 0)
            {
                b_targets[index]      = targets_tmp[op_perm[type][i]];
                b_operators[index]    = b_operators_tmp[op_perm[type][i]];
                b_target_ids[index++] = target_ids_tmp[op_perm[type][i]];
                number_in_groups++;
            }
            if (type == 1)
            {
                size_t source_index = op_perm[type][i];
                for (size_t j = 0; j < 2; ++j)
                {
                    b_targets[index]      = targets_tmp[source_index];
                    b_operators[index]    = b_operators_tmp[source_index];
                    b_target_ids[index++] = target_ids_tmp[source_index++];
                    number_in_groups++;
                }
            }

            if (type == 2)
            {
                size_t source_index = op_perm[type][i];
                for (size_t j = 0; j < 4; ++j)
                {
                    b_targets[index]      = targets_tmp[source_index];
                    b_operators[index]    = b_operators_tmp[source_index];
                    b_target_ids[index++] = target_ids_tmp[source_index++];
                    number_in_groups++;
                }
            }
            if (type == 3)
            {
                size_t source_index = op_perm[type][i];
                for (size_t j = 0; j < 8; ++j)
                {
                    b_targets[index]      = targets_tmp[source_index];
                    b_operators[index]    = b_operators_tmp[source_index];
                    b_target_ids[index++] = target_ids_tmp[source_index++];
                    number_in_groups++;
                }
            }
        }
    }
}

template <typename Real, typename arch>
void Box<Real, arch>::set_omega(CoeffMatrix *o)
{
    omega = o;
}

template <typename Real, typename arch>
void Box<Real, arch>::set_mu(CoeffMatrix *m)
{
    mu = m;
}

template <typename Real, typename arch>
void Box<Real, arch>::set_offsets(size_t local_index, size_t offset_id, Real3 periodic_shift /*= Real3(0.,0.,0.)*/)
{
    particle_offset_ids[local_index]      = offset_id;
    particle_periodic_shifts[local_index] = periodic_shift;
}

template <typename Real, typename arch>
void Box<Real, arch>::set_a_interaction_pairs(size_t local_index, CoeffMatrix *omega, M2M_Operator* global_a_ptr, size_t omega_id)
{
    a_operators[local_index]  = global_a_ptr;
    a_targets[local_index]    = omega;
    a_target_ids[local_index] = omega_id;
}

template <typename Real, typename arch>
void Box<Real, arch>::set_b_interaction_pairs(size_t local_index, CoeffMatrix *mu, M2L_Operator* global_b_ptr, size_t mu_id)
{
    b_operators[local_index]  = global_b_ptr;
    b_targets[local_index]    = mu;
    b_target_ids[local_index] = mu_id;
}

template <typename Real, typename arch>
void Box<Real, arch>::set_c_interaction_pairs(size_t local_index, CoeffMatrix *mu, L2L_Operator* global_c_ptr, size_t mu_id)
{
    c_operators[local_index]  = global_c_ptr;
    c_targets[local_index]    = mu;
    c_target_ids[local_index] = mu_id;
}

template class Box<double, Device<double> >;
template class Box<float, Device<float> >;

template class Box<double, Host<double> >;
template class Box<float, Host<float> >;

} //namespace end
