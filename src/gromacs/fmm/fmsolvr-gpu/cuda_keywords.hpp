#ifndef _BK_cuda_keywords_hpp_
#define _BK_cuda_keywords_hpp_
#include <cstdio>

#ifdef __CUDACC__
#define CUDA   __host__ __device__
#define DEVICE __device__
#else
#define CUDA
#define DEVICE
#endif  

// Define this to turn on error checking
#define CUDA_ERROR_CHECK

#ifdef __CUDACC__

#define CUDA_SAFE_CALL( err ) __cuda_safe_call( err, __FILE__, __LINE__ )
#define CUDA_CHECK_ERROR() __cuda_check_error( __FILE__, __LINE__ )
 
inline void __cuda_safe_call( cudaError err, const char *file, const int line )
{
  #ifdef CUDA_ERROR_CHECK
  if ( cudaSuccess != err )
  {
    fprintf( stderr, "cudaSafeCall() failed at %s:%i : %s\n",file, line, cudaGetErrorString( err ) );
    exit( -1 );
  }
  #endif
 
  return;
}
 
inline void __cuda_check_error( const char *file, const int line )
{
  #ifdef CUDA_ERROR_CHECK
  cudaError err = cudaGetLastError();
  if ( cudaSuccess != err )
  {
    fprintf( stderr, "cudaCheckError() failed at %s:%i : %s\n",file, line, cudaGetErrorString( err ) );
    exit( -1 );
  }
  #endif
  return;
}
#endif


#endif
