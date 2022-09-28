#ifndef _BK_fmm_hpp_
#define _BK_fmm_hpp_

#define P2P_STREAMS 4
#define STREAMS 8

#include <thread>
#include <array>

#include "xyzq.hpp"
#include "cuda_alloc.hpp"
#include "box.hpp"
#include "input_output.hpp"
#include "fmm_printer.hpp"
#include "global_functions.hpp"
#include "architecture.hpp"
#include "multipole2multipole.hpp"
#include "local2local.hpp"
#include "input_output.hpp"
#include "ioadapter.hpp"
#include "architecture.hpp"
#include "timer.hpp"
#include "latticeoperator.hpp"
#include "bit_manipulator.hpp"
#include "testdata.hpp"
#include "dipole_compensation.hpp"

namespace gmx_gpu_fmm{

class fmm_type{
public:
    typedef REAL Real;
};
/*!
 * \brief The fmm_algorithm handles the complete FMM computation. It implements all FMM stages.
 */
class fmm_algorithm : public fmm_type {

public:
    typedef fmm_type::Real Real;

    //! Class for handling the unsorted input particles, input box size and periodic conditions.
    //! It also provides the output structures for forces and potentials at particle positions
    typedef input_output<Real> io_type;
    //! Datatype for vector field of forces
    typedef typename io_type::field_vector_type           field_vector_type;
    //! Datatype for particles (4D vector field)
    typedef typename io_type::particle_vector_type        particle_vector_type;
    //! Datatype for potentilas (scalar field)
    typedef typename io_type::potential_vector_type       potential_vector_type;
    //! Collective structure for the result on Device (forces and potentials)
    typedef typename io_type::outputadapter_type          outputadapter_type;
    //! Collective structure for the result on th Host (forces and potentials)
    typedef typename io_type::host_outputadapter_type     host_outputadapter_type;
    //typedef typename io_type::inputadapter_type         inputadapter_type;
    //! Datatype for particles (4D vector field) on the Host
    typedef typename io_type::host_particle_vector_type   host_particle_vector_type;
    //! Unique pointer to the {@link input_output} datasructure
    std::unique_ptr<io_type> io_unique;
    //! Pointer to the {@link input_output} datasructure
    io_type* io;

    //! Basic 3D datatype
    typedef XYZ<Real>  Real3;
    //! Basic 4D datatype
    typedef XYZQ<Real> Real4;
    //! Basic dataype for basis vectors describing the simulation box
    typedef ABC<Real3> Real33;

    //! Architecture (Device double) abstraction class for memory allocations
    typedef Device<double> DeviceD;
    //! Architecture (Device real) abstraction class for memory allocations
    typedef Device<Real> Device;
    //! Architecture (Device real) abstraction class for memory allocations
    typedef Host<Real> Host;
    //! Architecture (Device real) abstraction class for memory allocations
    typedef DeviceOnly<Real> DeviceOnly;

    //! Device allocator type for 4D datastructure
    typedef typename Device::allocator::template rebind<Real4>::other alloc_real4;
    //! Helper vector for dipole compenstion
    std::vector<Real4, alloc_real4 > q0abc;
    //! Helper vector for dipole compenstion
    std::vector<Real4, alloc_real4 > fake_particles;

    //! Dipole compenstion Handler class type
    typedef typename fmsolvr::dipole_compensation<fmm_algorithm> DipoleCompensation;
    //! Pointer to Dipole compensation implementation object
    DipoleCompensation* Dipole_compensation;
    //! Helper host vector for dipole compenstion
    std::vector<Real4, cuda_host_allocator<Real4> > q0abc_host;
    //! Helper host vector for dipole compenstion
    std::vector<Real4, cuda_host_allocator<Real4> > fake_particles_host;
    //! Number of fake particles, helper for dipole compenstion
    size_t fake_particle_size;

    //! Basic structre for handlig multipoles on Device, triangular array
    typedef MultipoleCoefficientsUpper<Real, Device> CoeffMatrix;
    //! Basic structre for handlig multipoles on Device, triangular array (double precision)
    typedef MultipoleCoefficientsUpper<double, DeviceD> CoeffMatrixD;
    //! Memory pool for {@link CoeffMatrix} to use allocation more efficiently
    typedef MultipoleCoefficientsUpper_Memory<Real, Device> CoeffMatrixMemory;
    //! Basic structre for handlig multipoles on the host, triangular array
    typedef MultipoleCoefficientsUpper<Real, Host> CoeffMatrix_Host;
    //! Memory pool for {@link CoeffMatrix_Host} to use allocation more efficiently
    typedef MultipoleCoefficientsUpper_Memory<Real, Host> CoeffMatrixMemory_Host;
    //! Basic structre for handlig multipoles, triangular array (double precision) for opertors
    typedef MultipoleCoefficientsUpper<double, DeviceD > OperatorMatrixDouble;

    //! Basic structre for handlig M2M operator
    typedef M2M_Operator<CoeffMatrix, OperatorMatrixDouble> M2M_Operator;
    //! Basic structre for handlig M2L operator
    typedef M2L_Operator<CoeffMatrix, OperatorMatrixDouble> M2L_Operator;
    //! Basic structre for handlig L2L operator
    typedef L2L_Operator<CoeffMatrix, OperatorMatrixDouble> L2L_Operator;

    CoeffMatrixMemory* omega_memory;
    CoeffMatrixMemory* mu_memory;
    CoeffMatrix* omegav;
    CoeffMatrix* muv;

    //! Array of pointers to Multipoles stored as {@link CoeffMatrix}
    CoeffMatrix** omega;
    //! Array of pointers to Taylor moments stored as {@link CoeffMatrix}
    CoeffMatrix** mu;

    CoeffMatrixMemory_Host* omega_memory_host;
    CoeffMatrixMemory_Host* mu_memory_host;

    CoeffMatrix_Host* omegav_host;
    CoeffMatrix_Host* muv_host;

    //! Array of pointers to Multipoles o stored as {@link CoeffMatrix_Host}
    CoeffMatrix_Host** omega_host;
    //! Array of pointers to Taylor moments stored as {@link CoeffMatrix_Host}
    CoeffMatrix_Host** mu_host;

    //! Target Taylormoment for receiving all periodic contributions in case of non periodic execution
    CoeffMatrix* non_periodic_mu_dummy;
    //! M2L operator for non periodic execution
    M2L_Operator* non_periodic_m2l_dummy;

    //! Basic type for storing the multipoles and Taylormoments as SoA structure
    typedef MultipoleCoefficientsUpperSoA<Real,Device> CoeffMatrixSoA;
    //! Pointer to SoA structure for multipoles
    CoeffMatrixSoA *omegaSoA;
    //! Pointer to SoA structure for taylormoments
    CoeffMatrixSoA *muSoA;

    //! Memory size for particles occupying boxes at lowest level. Checked in init phase
    std::vector<size_t> initial_box_particle_mem_sizes;
    //! Memory for in box local offsets of direct neighbors particles
    size_t* box_particle_source_offsets;
    //! Memory for in box local offsets of direct neighbors particles
    size_t* box_particle_source_starts;
    //! Memory for in box local offsets of direct neighbors particles
    size_t* box_particle_source_ends;
    //! Memory for in box stored periodic shifts of direct neighbors
    Real3* box_particle_source_shifts;

    //! Box datatype
    typedef Box<Real,Device> Box;
    //! Array of of octree boxes
    Box* box;
    //! Memory for {@link box} M2M operators list
    M2M_Operator** box_a_operators;
    //! Memory for {@link box} M2M targets list
    CoeffMatrix** box_a_targets;
    size_t* box_a_target_ids;
    //! Memory for {@link box} M2L operators list
    M2L_Operator** box_b_operators;
    //! Memory for {@link box} M2L targets list
    CoeffMatrix** box_b_targets;
    size_t* box_b_target_ids;
    //! Memory for {@link box} L2L operators list
    L2L_Operator** box_c_operators;
    //! Memory for {@link box} L2L targets list
    CoeffMatrix** box_c_targets;
    size_t* box_c_target_ids;
    //! Memory for {@link box} M2L operator temp for resorting
    M2L_Operator** box_operators_tmp;
    //! Memory for {@link box} M2L targets temp for resorting
    CoeffMatrix** box_targets_tmp;
    size_t* box_target_ids_tmp;

    typedef typename CoeffMatrix::value_type value_type;
    //! Resorted pointers to SoA targets for M2M operator
    value_type*** box_a_targets_SoA;
    //! Resorted pointers to SoA targets for M2L operator
    value_type*** box_b_targets_SoA;
    //! Resorted pointers to SoA targets for L2L operator
    value_type*** box_c_targets_SoA;
    //! Memory for Taylormoments derivatives in x direction
    CoeffMatrixMemory* dmux_memory;
    //! Memory for Taylormoments derivatives in y direction
    CoeffMatrixMemory* dmuy_memory;
    //! Memory for Taylormoments derivatives in z direction
    CoeffMatrixMemory* dmuz_memory;
    //! Taylormoments derivatives in x direction
    CoeffMatrix **dmux;
    //! Taylormoments derivatives in y direction
    CoeffMatrix **dmuy;
    //! Taylormoments derivatives in z direction
    CoeffMatrix **dmuz;

    CoeffMatrix* dmuxv;
    CoeffMatrix* dmuyv;
    CoeffMatrix* dmuzv;

    //! Pointer of pointers to M2M operator pointers (depth, direction)
    M2L_Operator *** m2l_operators;
    //! Pointer of pointers to M2L operator pointers (depth, direction)
    M2M_Operator *** m2m_operators;
    //! Pointer of pointers to L2L operator pointers (depth, direction)
    L2L_Operator *** l2l_operators;
    //! Lattice operator
    CoeffMatrix* Lattice;
    //! Lattice operator in double precision
    CoeffMatrixD * LatticeD;
    //! scaling factor for lattice operator
    double lattice_scale;
    //! rescaling factor for lattice operator
    Real lattice_rescale;
    //! list of dummy boxes for each depth
    M2L_Operator** dummy_m2l_operators;

    //! Multipole order
    const size_t p;
    //! Separation criterion
    const ssize_t ws;
    //! Depth of the octree
    size_t depth;
    //! Open or periodic
    const bool open_boundary_conditions;
    //! Compensate dipole
    const bool dipole_compensation;
    //! dipole of the simulation box
    Real3 dipole_host;
    Real3* dipole;
    //! spaese fmm type
    const bool sparse_fmm;

    int p2p_version;

    //! {@link p} + 1
    size_t p1;
    //! {@link p1} * {@link p1}
    size_t p1xx2;
    //! ({@link p } + 1) * ({@link p } + 2) / 2
    size_t p1xp2_2;
    size_t pxp1_2;
    //! NUmber of all boxes in the octree
    size_t num_boxes_tree;
    //! {@link p } * ({@link p } + 1) / 2
    //! NUmber of all boxes at depth={@link depth}
    size_t num_boxes_lowest;
    //! NUmber of all boxes above depth={@link depth}
    size_t global_offset;
    //! Number of empty boxes (no particles)
    size_t empty_boxes;
    //! (2*{@link ws}+1)*2;
    ssize_t ws_dim;
    //! {@link ws_dim}*{@link ws_dim}
    ssize_t ws_dim2;
    //! {@link ws_dim}*{@link ws_dim}*{@link ws_dim}
    ssize_t ws_dim3;
    //! CUDA blocksize for P2M annd Forces kernel
    size_t bs;
    //! CUDA number of blocks for P2M annd Forces kernel
    size_t n_blocks;
    //! Number of all M2L possible M2L operators
    size_t num_of_all_ops;
    //! Number of all effective M2L operators due to sibling position assymetry
    size_t num_of_efective_ops;
    //! Global maximal number of particles residnig in all boxes on the lowest level
    size_t max_particles_in_box;

    typedef typename Device::allocator::template rebind<size_t>::other alloc_size_t;

    //! CUDA blocks map to box ids
    std::vector<size_t, alloc_size_t> block_map;
    //! CUDA cummulative offset emerging from mapping box particles to fixed CUDA block size
    std::vector<size_t, alloc_size_t> offset_map;
    //! Mapping of octree lexicogrpahical ids to linear {@link box} order
    std::vector<size_t, alloc_size_t> box_id_map;
    //! Original ids of resorted particles
    std::vector<size_t, alloc_size_t> orig_ids;
    //! Ordered ids of originallz sorted particles
    std::vector<size_t, alloc_size_t> fmm_ids;
    //! CUDA blocks map to box ids on the Host
    std::vector<size_t, cuda_host_allocator<size_t> > block_map_host;
    //! CUDA cummulative offset emerging from mapping box particles to fixed CUDA block size on the Host
    std::vector<size_t, cuda_host_allocator<size_t> > offset_map_host;

    typedef typename Device::allocator::template rebind<Real3>::other alloc_real3;
    //! Octree boxes centers for multipole expansion origins
    std::vector<Real3, alloc_real3> expansion_points;
    //! Instance of printer class for debugging and output
    fmm_printer<fmm_algorithm> printer;

    typedef uint32_t boxid_type;
    typedef Real4 particle_type;

    //! Array of ordered particles, distributed to octree boxes
    particle_vector_type ordered_particles;
    //! Array of ordered particles at host, distributed to octree boxes
    host_particle_vector_type ordered_particles_host;
    //! Array of particle offstes
    std::vector<size_t, alloc_size_t> box_particle_offset;
    //! Array of particle offstes at host
    std::vector<size_t, cuda_host_allocator<size_t> > box_particle_offset_host;

    //int2 * sparse_interactions_list;
    //int* sparse_interactions_list_size;

    typedef typename Device::allocator::template rebind<int>::other alloc_int;
    //! Number of particles per box for CUDA atomics
    std::vector<int, alloc_int> particles_per_box_int;
    //! Array of particle offstes as integer for CUDA atomics
    std::vector<int, alloc_int> box_particle_offset_int;
    //! Memory for exclusion map
    int*  exclusions_memory;
    //! Interactions to be excluded
    int** exclusions;
    //! Sizes of exclusion maps
    int*  exclusions_sizes;
    //! exclusions pairs
    int2* exclusion_pairs;
    //! exclusions pairs
    int exclusion_pairs_size;
    //! Pointer to bitset to M2L operator
    Bitit<32,5,Device>* bitset;
    //! Memory for bitset
    unsigned int* bitset_mem;
    //! Object of a time measuring class
    Exe_time exec_time;
    //! P2P box pairs
    int2 **  p2p_particles_box_pairs;
    //! P2P box pairs on the Host
    int2 **  p2p_particles_box_pairs_host;
    //! P2P box particle sizes / 8
    int*  p2p_particles_sizes_div8;
    //! P2P box particle sizes rounded to mod32 = 0
    int*  p2p_rounded_particles_sizes;
    //! P2P particles
    REAL4**  p2p_particles;
    //! P2P periodic shifts for particles
    REAL3**  p2p_particles_periodic_shifts;
    //! P2P periodic shifts for particles on the Host
    REAL3**  p2p_particles_periodic_shifts_host;
    //! P2P forces and potentials at {@link p2p_particles} positions
    REAL4**  p2p_results;
    //! P2P particles on the host
    Real4** p2p_particles_host;
    //! P2P particles results on the host
    Real4** p2p_results_host;
    //! CUDA streams for P2P
    std::array<cudaStream_t, P2P_STREAMS> streams;
    //! CUDA streams for Farfield
    std::array<cudaStream_t, STREAMS> priority_streams;
    //! CUDA events for synchronizing for P2P
    std::array<cudaEvent_t, P2P_STREAMS> events;
    //! CUDA events for synchronizing for Farfield
    std::array<cudaEvent_t, STREAMS> priority_events;
    //! CUDA event
    cudaEvent_t copy_particles_to_host_event;
    //! CUDA event
    cudaEvent_t ordered_particles_set_event;
    //!Energy chunks
    Real* Ec;
    //!lenard Jones parameter
    REAL2* lj_comb;


    //cudaStream_t streams[P2P_STREAMS];
    //cudaStream_t priority_streams[STREAMS];

    //cudaEvent_t events[P2P_STREAMS];
    //cudaEvent_t priority_events[STREAMS];
    //! CUDA Stream counter
    int current_stream;
    //! CUDA Stream counter
    int current_priority_stream;

    //! Constructor
    fmm_algorithm();

    /*!
     * \brief fmm_algorithm
     * \param p                  Multipole moment
     * \param ws                 Separation criterion
     * \param depth              Octree depth
     * \param open_boundary      Openboundary flag
     * \param compensate_dipole  Dipole compensation flag
     */
    fmm_algorithm(size_t p, ssize_t ws, size_t depth, bool open_boundary, bool compensate_dipole, bool sparse);

    /*!
     * \brief Allocates memory for the result pointer of the {@link io}
     */
    void alloc_result();
    /*!
     * \brief Initializes the {@link io} object
     * \param n_particles Number of particles
     */
    void init_io(size_t n_particles);
    /*!
     * \brief Scales the simulation box
     * \param abc Basis vectors describint the box geometry (only cubic possible)
     */
    void update_fmm_boxscale(Real33 &abc);

    ~fmm_algorithm();
    /*!
     * \brief Allocates memory and sets exclusions
     * \param excl Exclusion map
     */
    void alloc_exclusions(std::vector<std::vector<int> > &excl);
    /*!
     * \brief Allocates memory for exclusions
     */
    void alloc_exclusions();
    /*!
     * \brief Distributes particles on the host at the precomputing phase to estimate the partice distribution
     * \param abc                     Simulationbox basis vectors
     * \param ref_corner              Reference corner
     * \param input                   Particle data
     * \param fmm_depth               Depth of the octree
     * \param[out] mem_size_per_box   Memory size needed per each box
     * \param n                       Number of all particles
     */
    void distribute_particle_memory(
            const Real33 & abc,
            const Real3 & ref_corner,
            const Real4* input,
            int fmm_depth,
            size_t* mem_size_per_box,
            size_t n);

    /*!
     * \brief Allocate memory for particles
     */
    void alloc_box_particle_memory(float factor);
    /*!
     * \brief Compute dipole compensating particles
     */
    void compute_dipole_compensating_particles();

    /*!
     * \brief Copy input particles to the {@link io}
     * \tparam Inputdata type
     * \param  Particles
     */
    template <typename Inputdata>
    void copy_particles(Inputdata& particles);
    /*!
     * \brief Distribute particles strater
     */
    void distribute_particles(bool gmx_does_neighbor_search = true);
    /*!
     * \brief The actual implementation of the particles ditribution
     */
    void distribute_particles_impl(bool gmx_does_neighbor_search);

    void gmx_copy_particles_buffer(REAL *inputparticles, cudaEvent_t *gmx_h2d_ready_event);

    void gmx_copy_particles_buffer_impl(REAL *inputparticles, cudaEvent_t *gmx_h2d_ready_event);
    /*!
     * \brief Set expansion points at each box in the tree
     */
    void set_expansion_points();
    /*!
     * \brief  Computes result energy on device
     * \return Energy
     */
    Real energy_impl();
    /*!
     * \brief Computes scaled result energy
     * \param eps  Epsilon
     * \return     Energy
     */
    double energy(Real eps = 1.0);
    /*!
     * \brief Computes energy on originally oredered particles
     * \return  Energy
     */
    double energy_orig_order();
    /*!
     * \brief prins result for each particle in original input order
     * \return  void
     */
    void dump_result_in_orig_order();
    /*!
     * \brief Compute and displys energy on the Device
     */
    void energy_dump();
    /*!
     * \brief Actual implementation of energy comutation on the Device
     * \param result_ptr   Pointer to the result structre
     */
    void energy_dump_impl(outputadapter_type* result_ptr);
    /*!
     * \brief  Compute and displys forces on the Device
     * \param type   Type of forces to dump, 0-all, 1-x, 2-y, 3-z direction
     */
    void force_dump(int type);
    /*!
     * \brief Actual implementation of forces comutation on the Device
     * \param result_ptr   Pointer to the result structre
     * \param type         Type of forces to dump, 0-all, 1-x, 2-y, 3-z direction
     */
    void force_dump_impl(outputadapter_type* result_ptr, int type);

    /*!
     * \brief Computes l2 norm of the forces
     * \return l2 norm
     */
    Real3 force_l2_norm();
    /*!
     * \brief Computes l1 norm of the forces
     * \return l1 norm
     */
    Real3 force_l1_norm();
    /*!
     * \brief Instantiates and precomputes M2M operators
     */
    void set_m2m_operator();
    /*!
     * \brief Instantiates and precomputes M2L operators
     */
    void set_m2l_operator();
    /*!
     * \brief Instantiates and precomputes L2L operators
     */
    void set_l2l_operator();
    /*!
     * \brief Allocates all data on Device and Host
     */
    void alloc_data();
    /*!
     * \brief Actual implementation of data allocation on Device and Host
     */
    void alloc_data_impl();
    /*!
     * \brief Frees all memory
     */
    void free_data();
    /*!
     * \brief Precomputes all needed data to execute P2P
     */
    void prepare_p2p(int p2p_version = 1);
    /*!
     * \brief Actual implementation of the P2P preparation
     */
    void prepare_p2p_impl(int version);
    /*!
     * \brief Precomputes all data nedded to execute M2M translation
     */
    void prepare_m2m();
    /*!
     * \brief Actual implementation of the M2M preparation
     */
    void prepare_m2m_impl();
    /*!
     * \brief Precomputes all data nedded to execute M2L translation
     */
    void prepare_m2l();
    /*!
     * \brief Actual implementation of the M2M preparation
     */
    void prepare_m2l_impl();
    /*!
     * \brief Precomputes all data nedded to execute L2L translation
     */
    void prepare_l2l();
    /*!
     * \brief Actual implementation of the M2M preparation
     */
    void prepare_l2l_impl();
    /*!
     * \brief Precomputes all data nedded to execute P2M and FORCES calculations
     */
    void prepare_data();
    /*!
     * \brief Actual implementation of prepare data function
     */
    void prepare_data_impl();
    /*!
     * \brief Calculates nearfield forces and potentials
     * \param host_device   Splits calculation to Host and Device. 13 all interactions computed on device, 0 all interactions computed on host
     * \param version       0 - old version, 1 - new version
     */
    void p2p(bool calc_energies = true, int host_device = 13, int version = 1);
    /*!
     * \brief Actual implementation of the P2P calculation
     * \param result_ptr   Pointer to the result structure
     */
    void p2p_impl(outputadapter_type* result_ptr);
    /*!
     * \brief Implementaion of the P2P in open boundary conditons
     * \param result_ptr   Pointer to the result structure
     */
    void p2p_impl_open(outputadapter_type* result_ptr);
    /*!
     * \brief Implementaion of the P2P in periodic boundary conditons
     * \param result_ptr   Pointer to the result structure
     * \param host_device  Host device flag
     * \param version      P2P algorithm version (0,1)
     */
    void p2p_impl_periodic(outputadapter_type* result_ptr, bool calc_energies, int host_device, int version);
    /*!
     * \brief Calculates nearfield forces and potentials on the Host
     * \param host_device  Host device flag
     * \param type         P2P algorithm version (0,1)
     */
    void p2p_host(int host_device, int type);
    /*!
     * \brief Implementaion of the P2P in open boundary conditons
     * \param host_device  Host device flag
     * \param type         P2P algorithm version (0,1)
     */
    void p2p_host_impl(int host_device, int type);
    /*!
     * \brief Executes the P2M stage of the FMM
     */
    void p2m();
    /*!
     * \brief Actual P2M implementation
     */
    void p2m_impl();
    /*!
     * \brief Executes the M2M stage of the FMM
     */
    void m2m();
    /*!
     * \brief Actual M2M implementation
     */
    void m2m_impl();
    /*!
     * \brief Executes the M2L stage of the FMM
     */
    void m2l();
    /*!
     * \brief Splits the M2L between Host and Device
     * \param box_chunks
     */
    void m2l_impl_gpu_cpu(std::vector<size_t>& box_chunks);

    //! Global box chunking for distribution between CPU and GPU
    std::vector<size_t> m2l_impl_gpu_cpu_prep();
    /*!
     * \brief Actual M2L implementation on Host
     */
    void m2l_impl_cpu();
    /*!
     * \brief Actual M2L implementation
     */
    void m2l_impl();
    /*!
     * \brief Executes the lattice computation on the GPU
     */
    void lattice();
    /*!
     * \brief Actual lattice implementation
     */
    void lattice_impl();
    /*!
     * \brief Executes the L2L stage of the FMM
     */
    void l2l();
    /*!
     * \brief Actual L2L implementation
     */
    void l2l_impl();
    /*!
     * \brief Executes the forces and exclusions computation
     */
    void forces(Real one4pieps0 = 1.0);
    /*!
     * \brief Actual implementation of forces and exclusions computation
     * \param result_ptr  Pointer to the result structure
     */
    void forces_impl(outputadapter_type* result_ptr, Real one4pieps0);
    /*!
     * \brief Copy forces back to the GROMACS forces container
     */
    void get_forces_in_orig_order(float3* gmx_forces, cudaStream_t *gmx_sync_stream);
    /*!
     * \brief Copy forces back to the GROMACS forces container
     */
    void get_forces_in_orig_order_impl(float3* gmx_forces, cudaStream_t *gmx_sync_stream);
    /*!
     * \brief Adds results from host and device if both used
     */
    void gather_results();
    /*!
     * \brief Calcualtes the L2 relative error norm of the taylormoments at the lowest level comared to given reference
     * \param sol   Array of pointers to FMM taylormoments
     * \param ref   Array of pointers to reference taylormoments
     * \return      relatice error in real and imag part
     */
    std::pair<double,double> get_error(CoeffMatrix** sol,CoeffMatrix** ref);

    void test_lj();
};

}//namespace end

#endif
