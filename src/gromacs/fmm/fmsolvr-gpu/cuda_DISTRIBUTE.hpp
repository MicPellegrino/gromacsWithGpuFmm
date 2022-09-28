#ifndef _BK_FMM_cuda_DISTRIBUTE_hpp
#define _BK_FMM_cuda_DISTRIBUTE_hpp

namespace gmx_gpu_fmm{


/*!
 * \brief Computes lexicographical 1D index based on 3D box position in the cubic octree. Device only function.
 * \param x       Position in x direction
 * \param y       Position in y direction
 * \param z       Position in z direction
 * \param depth   Depth of the octree
 * \return        size_t
 */
__device__
extern size_t __make_boxid(size_t x, size_t y, size_t z, unsigned depth)
{
    return (z << (depth * 2)) + (y << depth) + x;
}

/*!
 * \brief  Sets the particles and its original id's
 * \tparam Real4       4D data structure representing particle (x,y,z,q).
 * \param  particles   Array of particles.
 * \param  orig_ids    Array of original global ids of the particles
 * \param  n           Number of particles
 */
template <typename Real4>
__global__
void __flush_particles(Real4 *particles, size_t *orig_ids, size_t n)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i<n)
    {
        orig_ids[i] = 0;
        particles[i] = Real4(0.,0.,0.,0.);
    }
}

/*!
 * \brief Sets particle local particle index, particle_box_offset, box active flag, particle per box to zero.
 * \param box                    Pointer to a box.
 * \param parts_per_box_int      Number of particle in this box.
 * \param box_part_offset_int    Global offset of particles sorted into this box.
 * \param num_boxes_tree         Number of all boxes in the tree.
 * \param num_boxes_lowest       Number of boxes at the lowest level.
 */
template <typename Box>
__global__
void __zero_box_index(Box *box, int* parts_per_box_int, int* box_part_offset_int, size_t num_boxes_tree, size_t num_boxes_lowest)
{
    size_t box_id = blockIdx.x * blockDim.x + threadIdx.x;

    //in case of depth == 0
    if(box_id == 0)
    {
        box_part_offset_int[1]  = 0;
    }

    if(box_id >= num_boxes_tree)
        return;

    box[box_id].ptcl_index = 0;

    if(box_id >= num_boxes_tree - num_boxes_lowest)
    {
        box[box_id].active = 0;
    }

    if(box_id < num_boxes_lowest)
    {
        parts_per_box_int[box_id]  = 0;
        box_part_offset_int[box_id]  = 0;
    }

    if(box_id == num_boxes_lowest)
    {
        box_part_offset_int[box_id]  = 0;
    }
}

template <typename REAL3, typename Real4, typename Real>
__global__
void
__gmx2fmm_particle_copy(Real  *gmx_particles,
                        Real4 *fmm_particles,
                        Real   box_scale,
                        int    num_particles)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i >= num_particles)
        return;

    Real4  fmmp;
    Real   q         = *(reinterpret_cast<Real*>(fmm_particles + i) + 3);
    REAL3* gmxp_ptr  = reinterpret_cast<REAL3*>(gmx_particles + i*3);
    REAL3  gmxp      = *gmxp_ptr;
    fmmp.x           = gmxp.x * box_scale;
    fmmp.y           = gmxp.y * box_scale;
    fmmp.z           = gmxp.z * box_scale;
    fmmp.q           = q;
    fmm_particles[i] = fmmp;

    //printf("%d %e %e %e\n",i,  gmxp.x,  gmxp.y, gmxp.z);//, gmx_particles[i].w);
    //printf("%d %e %e %e %e\n",i,  fmm_particles[i].x,  fmm_particles[i].y, fmm_particles[i].z, fmm_particles[i].q);

}

/*!
 * \brief Distributes the the gloal particl indices to the boxes at the lowest level and computes the number of particles in the box.
 * Detemines the box id according to the (x,y,z) coordinates of the input particle.
 * \tparam Real33             3D basis vector type for simulations box.
 * \tparam Real3              Vector type for position.
 * \tparam Real4              Particle type for position and charge.
 * \tparam Box                Octree box type
 * \tparam Real               Underlying scalar datatype.
 * \param convert_to_abc      Standard basis.
 * \param ref_corner          Reference corner of the simulation box.
 * \param input               Array of particles.
 * \param parts_per_box_int   Global array string number of particles per box.
 * \param fmm_depth           Depth of the octree.
 * \param depth_offset        Number of boxes at boxes above depth fmm_depth.
 * \param box_id_map          Map of the octree lexicographical order ids to box pointers.
 * \param box                 Current box for computation.
 * \param scale               Box scale for computing relatice positions.
 * \param particles_tile      Number of particles to process => global n/number of streams.
 * \param particles_offset    Particle offset if run concurrenntly.
 */
template <typename Real33, typename Real3, typename Real4, typename Box, typename Real>
__global__
void
__distribute_particles(
        const Real33 abc,
        const Real33 convert_to_abc,
        const Real3 ref_corner,
        Real4* input,
        int* parts_per_box_int,
        int fmm_depth,
        size_t depth_offset,
        size_t* box_id_map,
        Box* box,
        Real scale,
        size_t particles_tile,
        size_t particles_offset)
{

    const size_t i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i < particles_tile)
    {
        typedef uint32_t boxid_type;

        boxid_type box_x;
        boxid_type box_y;
        boxid_type box_z;
        boxid_type boxid;
        size_t index = i + particles_offset;
        Real4 particle = input[index];
        Real3 position = Real3(particle);

        if (position.x >= abc.a.x)
        {
            position.x -= abc.a.x;
        }
        if(position.x < 0.0)
        {
            position.x += abc.a.x;
        }
        if (position.y >= abc.b.y)
        {
            position.y -= abc.b.y;
        }
        if(position.y < 0.0)
        {
            position.y += abc.b.y;
        }
        if (position.z >= abc.c.z)
        {
            position.z -= abc.c.z;
        }
        if(position.z < 0.0 )
        {
            position.z += abc.c.z;
        }

        Real3 normalized = convert_to_abc * ( position - ref_corner );
        
        box_x = normalized.x * scale;
        box_y = normalized.y * scale;
        box_z = normalized.z * scale;

        boxid = __make_boxid(box_x, box_y, box_z, fmm_depth);

        box[box_id_map[depth_offset + boxid]].set_orig_index(index);
        atomicAdd(&parts_per_box_int[boxid],1);

        particle.x   = position.x;
        particle.y   = position.y;
        particle.z   = position.z;
        input[index] = particle;

    }
}

__global__
void
__compute_offset(int* parts_per_box_int,
                 int* box_part_offset_int,
                 size_t num_boxes_lowest)
{

    const size_t target_id = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t source_id = blockIdx.y * gridDim.z + blockIdx.z;

    if(target_id > num_boxes_lowest)
        return;
    if(source_id >= num_boxes_lowest)
        return;

    if(target_id == 0)
    {
        box_part_offset_int[target_id] = 0;
    }
    else if(source_id < target_id)
    {
        atomicAdd(&box_part_offset_int[target_id], parts_per_box_int[source_id]);
    }
}

__global__
void
__recast_offset( size_t* box_part_offset,
                 int* box_part_offset_int,
                 size_t num_boxes_lowest)
{

    const size_t id = blockIdx.x * blockDim.x + threadIdx.x;

    if(id >= num_boxes_lowest+1)
        return;

    box_part_offset[id] = (size_t)box_part_offset_int[id];

}

}//namespace end



#endif
