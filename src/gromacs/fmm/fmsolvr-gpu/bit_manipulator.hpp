#ifndef _BK_bit_manipulator_hpp
#define _BK_bit_manipulator_hpp

#include <bitset>
#include "managed.hpp"
#include "architecture.hpp"
#include "cuda_keywords.hpp"

namespace gmx_gpu_fmm{

/*! \brief Implements the bitset on Device that is used by the M2L operator for stroing the signs of its values.
 *  \tparam size_of_t Number of bits of the underlying datastructre.
 *  \tparam pow_of_2  The exponent pow_of_2 such that \f$2^{pow\_of\_2} = size\_of\_t\f$.
 *  \tparam arch      Underlying architecture [Host or Device].
*/
template <size_t size_of_t, size_t pow_of_2, typename arch>
class Bitit : public Managed<arch>{

public:

    //! Value type the bitset is based on. Must match the {@link size_of_t} in size of bytes.
    typedef unsigned int value_type;

    //! Stores value of 1. Helper for operations on the bit level.
    value_type one;
    //! Stores value of {@link size_of_t} - 1. Helper for operations on bit level.
    value_type size_of_t_minus_1;
    //! Memory for the bitset.
    value_type* bits;
    //! Size of the bitset == length of {@link bits}.
    value_type size;

    //! Empty constructor.
    Bitit();

    /*! \brief Initializes an empty Bitset.
        \param bits_ptr    Pointer to the bitset memory.
        \param sizeofarray Number of needed value_types to store n bits.
    */
    void init(value_type* bits_ptr, value_type sizeofarray );

    //! Initializes the Bitset with random values.
    //! Only for testing purposes.
    void set_random(value_type seed);

    //! Prints the bitset
    void dump();

    /*! \brief Computes an offset of elements of {@link bits} to get bit at position \p pos.
        \param pos Position of the bit in the bitset.
    */
    CUDA
    value_type offset(value_type pos);

    /*! \brief Computes bitmask to extract the bit at position \p pos
        \param pos Position of the bit in the bitset.
    */
    CUDA
    value_type mask(value_type pos);

    /*! \brief Sets the at position \p pos to 1. Position counting in array like manner.
        \param pos Position of the bit in the bitset.
    */
    void set(value_type pos);

    /*! \brief Sets the at position \p pos to 0. Position counting in array like manner.
        \param pos Position of the bit in the bitset.
    */
    void reset(value_type pos);

    /*! \brief Sets the at position \p pos to 1. Position counting in array like manner. Device function only.
        \param pos Position of the bit in the bitset.
    */
    DEVICE
    void cuset(value_type pos);

    /*! \brief Sets the at position \p pos to 0. Position counting in array like manner. Device function only.
        \param pos Position of the bit in the bitset.
    */
    DEVICE
    void cureset(value_type pos);

    /*! \brief  Returns the value of the bit at position \p pos. Position counting in array like manner. Device function only.
        \param  pos Position of the bit in the bitset.
        \return Boolean
    */
    CUDA
    bool operator ()(value_type pos);

    /*! \brief  Returns the value of the bit at position \p pos. Position counting in array like manner. Device function only.
        \param  pos Position of the bit in the bitset.
        \return Boolean
    */
    CUDA
    bool get(value_type pos);

    /*! \brief  Returns the element that stores the bit at positoin \p pos.
        \param  pos Position of the bit in the bitset.
        \return Reference to value_type.
    */
    CUDA
    value_type &get_item(value_type pos);
};

}//namespace end

#endif
