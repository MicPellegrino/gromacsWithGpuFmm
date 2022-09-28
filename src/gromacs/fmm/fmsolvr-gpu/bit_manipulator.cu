#include "bit_manipulator.hpp"

namespace gmx_gpu_fmm{

template <size_t size_of_t, size_t pow_of_2, typename arch>
Bitit< size_of_t,  pow_of_2,  arch>::Bitit() : size_of_t_minus_1(size_of_t-1),one(1){}

template <size_t size_of_t, size_t pow_of_2, typename arch>
void Bitit< size_of_t,  pow_of_2,  arch>::init(value_type* bits_ptr, value_type sizeofarray )
{
    bits = bits_ptr;
    size = sizeofarray;
    for (size_t i = 0; i < size; ++i)
    {
        bits[i] = 0;
    }
}

template <size_t size_of_t, size_t pow_of_2, typename arch>
void Bitit< size_of_t,  pow_of_2,  arch>::set_random(value_type seed)
{
    value_type s = seed;
    for(size_t i = 0; i < size; i++)
    {
        s *=(i+1);
        s ^=s <<13;
        s ^=s >>7;
        s ^=s <<17;
        bits[i] = s;
    }
}

template <size_t size_of_t, size_t pow_of_2, typename arch>
void Bitit< size_of_t,  pow_of_2,  arch>::dump()
{
    //printf("bits address %p\n",bits);
    //printf("size %d\n",size);
    for(size_t i = 0; i < size; i++)
    {
        if(i<10)
            printf("%lu  ",i);
        else
            printf("%lu ",i);
        std::bitset<size_of_t> bitset(bits[i]);
        for(size_t i = size_of_t; i>0; i--)
        {
                std::cout<<bitset[i-1];
        }
        printf("\n");
    }
    std::cout<<std::endl;
}

template <size_t size_of_t, size_t pow_of_2, typename arch>
CUDA
typename Bitit< size_of_t,  pow_of_2,  arch>::value_type  Bitit< size_of_t,  pow_of_2,  arch>::offset(value_type pos){
    return ( pos + (( pos >> size_of_t_minus_1) & (one<<pow_of_2 + ~0))) >> pow_of_2;
}

template <size_t size_of_t, size_t pow_of_2, typename arch>
CUDA
typename Bitit< size_of_t,  pow_of_2,  arch>::value_type  Bitit< size_of_t,  pow_of_2,  arch>::mask(value_type pos){
    return one<<(size_of_t_minus_1 - pos&size_of_t_minus_1);
}

//counting in array manner
//usual way
template <size_t size_of_t, size_t pow_of_2, typename arch>
void Bitit< size_of_t,  pow_of_2,  arch>::set(value_type pos)
{
    value_type off = offset(pos);
    bits[off] |= mask(pos);
}

template <size_t size_of_t, size_t pow_of_2, typename arch>
void Bitit< size_of_t,  pow_of_2,  arch>::reset(value_type pos)
{
    value_type off = offset(pos);
    bits[off] &= ~mask(pos);
}

template <size_t size_of_t, size_t pow_of_2, typename arch>
DEVICE
void Bitit< size_of_t,  pow_of_2,  arch>::cuset(value_type pos)
{
    value_type off = offset(pos);
    atomicOr(&bits[off],mask(pos));
}

template <size_t size_of_t, size_t pow_of_2, typename arch>
DEVICE
void Bitit< size_of_t,  pow_of_2,  arch>::cureset(value_type pos)
{
    value_type off = offset(pos);
    atomicAnd(&bits[off],~mask(pos));
}

template <size_t size_of_t, size_t pow_of_2, typename arch>
CUDA
bool Bitit< size_of_t,  pow_of_2,  arch>::operator ()(value_type pos)
{
    value_type off = offset(pos);
    return bits[off] & mask(pos);
}

template <size_t size_of_t, size_t pow_of_2, typename arch>
CUDA
bool Bitit< size_of_t,  pow_of_2,  arch>::get(value_type pos)
{
    value_type off = offset(pos);
    return bits[off] & mask(pos);
}

template <size_t size_of_t, size_t pow_of_2, typename arch>
CUDA
typename Bitit< size_of_t,  pow_of_2,  arch>::value_type& Bitit< size_of_t,  pow_of_2,  arch>:: get_item(value_type pos)
{
    value_type off = offset(pos);
    return bits[off];
}

template class Bitit<32,5, Device<float> >;
template class Bitit<32,5, Device<double> >;

}//namespace end
