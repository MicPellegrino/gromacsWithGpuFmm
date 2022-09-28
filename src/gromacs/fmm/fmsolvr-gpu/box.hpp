#ifndef BOX_HPP
#define BOX_HPP

#include <vector>
#include "managed.hpp"
#include "multipole.hpp"
#include "architecture.hpp"
#include "xyz.hpp"
#include "global_functions.hpp"
#include "multipole2multipole.hpp"
#include "multipole2local.hpp"
#include "local2local.hpp"

namespace gmx_gpu_fmm{

/*! \brief FMM tree Box class.
 * \tparam Real Basic scalar datatype
 * \tparam arch Underlying architecture
 * The class precomputes required data for operators at the octree box level.
 * It strores all information to perform valid far field and near field operations.
*/
template <typename Real, typename arch>
class Box : public Managed<arch>
{
public:
    //! The original ids of particles residing in this box.
    size_t* orig_ptcl_ids;

    //! The Multipole/Taylormoment data structre double type
    typedef MultipoleCoefficientsUpper<double,Device<double> > CoeffMatrixD;
    //! The Multipole/Taylormoment data structre current data type
    typedef MultipoleCoefficientsUpper<Real,arch> CoeffMatrix;
    //! The M2L operator data type
    typedef M2L_Operator<CoeffMatrix, CoeffMatrixD> M2L_Operator;
    //! The M2M operator data type
    typedef M2M_Operator<CoeffMatrix, CoeffMatrixD> M2M_Operator;
    //! The L2L operator data type
    typedef L2L_Operator<CoeffMatrix, CoeffMatrixD> L2L_Operator;

    typedef typename CoeffMatrix::complex_type complex_type;

    //! position vector data type
    //! \tparam Real Basic scalar datatype
    typedef XYZ<Real> Real3;

    //! Pointer to Multipole moment of the current Box
    CoeffMatrix* omega;
    //! Pointer to Taylor moment of the current Box
    CoeffMatrix* mu;
    //! FMM Octree box Id
    size_t id;
    //! FMM Octree depth
    size_t depth;
    //! FMM Octree depth offset = number of all boxes above the current deepst level
    size_t depth_offset;
    //! local index of particles = number of particles
    int ptcl_index;
    //! Set to zero if no particles in the box
    int active;
    //! An array of periodic shifts of all 26 direct neighbors
    //! If neighbor box not periodical, the value is (0,0,0)
    Real3* particle_periodic_shifts;
    //! The offsets of particles corresponding to the nearest boxes as in particle_periodic_shifts
    size_t* particle_offset_ids;
    //! The temporary array of pointers to M2L operators for shuffling only
    M2L_Operator** b_operators_tmp;
    //! The temporary array of targets of the M2L operation (mu's) for shuffling only
    CoeffMatrix** targets_tmp;

    size_t* target_ids_tmp;

    //! The array of pointers to all valid M2M operators
    M2M_Operator** a_operators;
    //! The array of pointers to all valid M2M targets (parent boxes). Preserves the order of a_operators
    CoeffMatrix** a_targets;
    size_t* a_target_ids;

    //! The array of pointers to all valid M2L operators
    M2L_Operator** b_operators;
    //! The array of pointers to all valid M2L targets. Preserves the order of b_operators
    CoeffMatrix** b_targets;
    size_t* b_target_ids;

    //! The array of pointers to all valid L2L operators
    L2L_Operator** c_operators;
    //! The array of pointers to all valid L2L targets (child boxes). Preserves the order of c_operators
    CoeffMatrix** c_targets;
    size_t* c_target_ids;

    //! Helper structure for resorting M2L operators
    std::vector<std::vector<size_t> > op_perm;

    //! Constructor
    Box();

    /*!
     * \brief Device function for setting the original global index of the current particle that is to be added.
     * \param orig_index Global index of the input particle.
     */
    DEVICE
    void set_orig_index(size_t orig_index);

    /*!
     * \brief Allocates memory for global indices indices.
     * \param n. Number of global indices == number of particles residing in this Box.
     */
    void alloc_mem(size_t n);

    //! Frees memory for all original indizes
    void free_mem();

    /*!
     * \brief Initializes an empty box.
     * \param id_ Octree box id.
     * \param d   Depth of the octree.
     */
    void init(size_t id_, size_t d);

    /*!
     * \brief Remaps pointer {@link particle_periodic_shifts} and {@link particle_offset_ids} to preallocated memory.
     * \param mem_index  Global preallocated position for {@link particle_offset_ids} and {@link particle_periodic_shifts}.
     * \param offset_mem Pointer to memory for {@link particle_offset_ids}.
     * \param shifts_mem Pointer to memory for {@link particle_periodic_shifts}.
     */
    void set_offset_mem(size_t mem_index, size_t* offset_mem, Real3* shifts_mem);

    /*!
     * \brief Remaps to preallocated memory for {@link a_operators} and {@link a_targets}.
     * \param op_mem_index  Global preallocated position for {@link a_operators} and {@link a_targets}.
     * \param o_mem         Pointer to memory for {@link a_operators}.
     * \param t_mem         Pointer to memory for {@link a_targets}.
     */
    void set_a_mem(size_t op_mem_index, M2M_Operator** o_mem, CoeffMatrix** t_mem, size_t* id_mem);

    /*!
     * \brief Remaps to preallocated memory for {@link b_operators} and {@link b_targets}.
     * \param op_mem_index  Global preallocated position for {@link b_operators} and {@link b_targets}.
     * \param o_mem         Pointer to memory for {@link b_operators}.
     * \param t_mem         Pointer to memory for {@link b_targets}.
     */
    void set_b_mem(size_t op_mem_index, M2L_Operator** o_mem, CoeffMatrix** t_mem, size_t* id_mem);

    /*!
     * \brief Remaps to preallocated memory for {@link c_operators} and {@link c_targets}.
     * \param op_mem_index  Global preallocated position for {@link c_operators} and {@link c_targets}.
     * \param o_mem         Pointer to memory for {@link c_operators}.
     * \param t_mem         Pointer to memory for {@link c_targets}.
     */
    void set_c_mem(size_t op_mem_index, L2L_Operator** o_mem, CoeffMatrix** t_mem, size_t* id_mem);

    /*!
    * \brief Remaps to external memory for {@link op_perm}
    * \param Global permutations memory pointer for {@link op_perm}
    */
    void set_permutations(std::vector<std::vector<size_t> >&permutations);

    /*!
     * \brief Permutes the M2L operator and target pointers for apprioprate Device access.
     * \param b_operators_temp      Operator temp memory for permutation calculations.
     * \param targets_temp          Target temp memory for permutation calculations.
     * \param num_of_efective_ops   Number of M2L operators to permute (ws=1 => 189).
     */
    void permute_ops(M2L_Operator** b_operators_temp, CoeffMatrix **targets_temp, size_t* target_ids_temp, size_t num_of_efective_ops);

    /*!
     * \brief Setter for {@link omega}
     * \param o  Pointer to CoeffMatrix
     */
    void set_omega(CoeffMatrix *o);

    /*!
     * \brief Setter for {@link mu}
     * \param o  Pointer to CoeffMatrix
     */
    void set_mu(CoeffMatrix *m);

    /*!
     * \brief Sets the values for {@link particle_periodic_shifts} and {@link particle_offset_ids}.
     * \param local_index     Local index of the neighbor box.
     * \param offset_id       Global octree Id of the neighbr box
     * \param periodic_shift  Value of the periodic shift of the neighbor box.
     */
    void set_offsets(size_t local_index, size_t offset_id, Real3 periodic_shift = Real3(0.,0.,0.));

    /*!
     * \brief Sets valid interaction lists for the M2M operator.
     * \param local_index   Local index of the interaction pair.
     * \param omega         The pointer to the interacting omega (parent of this box).
     * \param global_a_ptr  The pointer to the valid operator for this M2M translation.
     */
    void set_a_interaction_pairs(size_t local_index, CoeffMatrix *omega, M2M_Operator* global_a_ptr, size_t omega_id);

    /*!
     * \brief Sets valid interaction lists for the M2L operator.
     * \param local_index   Local index of the interaction pair.
     * \param mu            The pointer to the interacting mu (child of the nearest neighbor of this box parent).
     * \param global_a_ptr  The pointer to the valid operator for this M2L transformation.
     */
    void set_b_interaction_pairs(size_t local_index, CoeffMatrix *mu, M2L_Operator* global_b_ptr, size_t mu_id);

    /*!
     * \brief Sets valid interaction lists for the L2L operator.
     * \param local_index   Local index of the interaction pair.
     * \param mu            The pointer to the interacting mu (child of this box).
     * \param global_a_ptr  The pointer to the valid operator for this M2L transformation.
     */
    void set_c_interaction_pairs(size_t local_index, CoeffMatrix *mu, L2L_Operator* global_c_ptr, size_t mu_id);

};

}//namespace end


#endif // BOX_HPP
