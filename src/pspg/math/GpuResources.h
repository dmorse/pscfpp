#ifndef GPU_RESOURCES_H
#define GPU_RESOURCES_H

#include <cufft.h>

extern int THREADS_PER_BLOCK;
extern int NUMBER_OF_BLOCKS;
//#define SINGLE_PRECISION
#define DOUBLE_PRECISION


#ifdef SINGLE_PRECISION
typedef cufftReal cudaReal;
typedef cufftComplex cudaComplex;
typedef cufftReal hostReal;
typedef cufftComplex hostComplex;
#else
#ifdef DOUBLE_PRECISION
typedef cufftDoubleReal cudaReal;
typedef cufftDoubleComplex cudaComplex;
typedef cufftDoubleReal hostReal;
typedef cufftDoubleComplex hostComplex;
#endif
#endif


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__);}
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
  if (code != cudaSuccess)
    {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if(abort) exit(code);
    }
}

static __global__ void helmholtzHelper(cudaReal* result, const cudaReal* composition,
   const cudaReal* pressure, float chi, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] = composition[i] * composition[i] / chi - pressure[i];
   }
}

static __global__ void reformField(cudaReal* Wa, cudaReal* Wb,
   const cudaReal* pressureF,const cudaReal* compositionF, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      Wa[i] = pressureF[i] + compositionF[i];
      Wb[i] = pressureF[i] - compositionF[i];
   }
}

static __global__ void mcftsStepHelper(cudaReal* result, const cudaReal* A2, const cudaReal* sqrtSq, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
          result[i] += A2[i] * sqrtSq[i];
   }
}
static __global__ void mcftsScale(cudaReal* result, cudaReal scale, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
          result[i] = result[i] * 2 * scale - scale;
   }
}

static __global__ void pointWiseAdd(cudaReal* result, const cudaReal* rhs, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] += rhs[i];
   }
}

static __global__ void subtractUniform(cudaReal* result, cudaReal rhs, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] -= rhs;
   }
}

static __global__ void addUniform(cudaReal* result, cudaReal rhs, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] += rhs;
   }
}

static __global__ void pointWiseSubtract(cudaReal* result, const cudaReal* rhs, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] -= rhs[i];
   }
}

static __global__ void pointWiseSubtractFloat(cudaReal* result, const float rhs, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] -= rhs;
   }   
}

static __global__ void pointWiseBinarySubtract(const cudaReal* a, const cudaReal* b, cudaReal* result, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] = a[i] - b[i];
   }
}

static __global__ void pointWiseBinaryAdd(const cudaReal* a, const cudaReal* b, cudaReal* result, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] = a[i] + b[i];
   }
}

static __global__ void pointWiseBinaryMultiply(const cudaReal* a, const cudaReal* b, cudaReal* result, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] = a[i] * b[i];
   }
}

static __global__ void pointWiseAddScale(cudaReal* result, const cudaReal* rhs, double scale, int size)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] += scale * rhs[i];
   }
}

static __global__ void assignUniformReal(cudaReal* result, cudaReal uniform, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < size; i += nThreads) {
      result[i] = uniform;
   }
}

static __global__ void assignReal(cudaReal* result, const cudaReal* rhs, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < size; i += nThreads) {
      result[i] = rhs[i];
   }
}

static __global__ void inPlacePointwiseMul(cudaReal* a, const cudaReal* b, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < size; i += nThreads) {
      a[i] *= b[i];
   }
}

static __global__ void reductionSum(cudaReal* sum, const cudaReal* in, int size)
{
   // number of blocks cut in two to avoid inactive initial threads
   int tid = threadIdx.x;
   int bid = blockIdx.x;
   int idx = bid * (blockDim.x*2) + tid;
   
   if (idx < size) {
      // Shared memory holding area
      extern __shared__ cudaReal sdata[];

      // Copy data into shared, wait for all threads to finish
      sdata[tid] = in[idx] + in[idx+ blockDim.x];
      
      __syncthreads();

      // Make comparisons across the block of data, each thread handling 
      // one comparison across two data points with strided indices before 
      // syncing with each other and then making further comparisons.
      for (int stride = blockDim.x/2; stride > 0; stride/=2) {
         if (tid < stride) {
            sdata[tid] += sdata[tid+stride];
         }
         __syncthreads();
      }

      // Store the output of the threads in this block
      if (threadIdx.x == 0) {
         sum[bid] = sdata[0];
      }
   }
}

static __global__ void reductionMaxAbs(cudaReal* max, const cudaReal* in, int size) 
{
   // number of blocks cut in two to avoid inactive initial threads
   int tid = threadIdx.x;
   int bid = blockIdx.x;
   int idx = bid * (blockDim.x*2) + tid;
   
   
   volatile extern __shared__ cudaReal sdata[];
   // global mem load combined with first operation. load with fabs.
   cudaReal in0 = fabs(in[idx]);
   cudaReal in1 = fabs(in[idx + blockDim.x]);
   sdata[tid] = (in0 > in1) ? in0 : in1;
   
   __syncthreads();
   
   // reduction
   // data and block dimensions need to be a power of two
   for (int stride = blockDim.x / 2; stride > 0; stride /= 2) {
      if (tid < stride) {
         if (sdata[tid+stride] > sdata[tid]) {
               sdata[tid] = sdata[tid+stride];
            }
      }
      __syncthreads();
   }
      
   // one thread for each block stores that block's results in global memory
   if (tid == 0)
      max[bid] = sdata[0];
}

static __global__ void reductionMax(cudaReal* max, const cudaReal* in, int size)
{
   // number of blocks cut in two to avoid inactive initial threads
   int tid = threadIdx.x;
   int bid = blockIdx.x;
   int idx = bid * (blockDim.x*2) + tid;
   
   if (idx < size) {
      // Shared memory holding area
      extern __shared__ cudaReal sdata[];

      // Copy data into shared, wait for all threads to finish
      sdata[tid] = (in[idx] > in[idx+blockDim.x] ) ? in[idx] : in[idx + blockDim.x];
      
      __syncthreads();

      // Make comparisons across the block of data, each thread handling 
      // one comparison across two data points with strided indices before 
      // syncing with each other and then making further comparisons.
      for (int stride = blockDim.x/2; stride > 0; stride/=2) {
         if (tid < stride) {
            if (sdata[tid+stride] > sdata[tid]) {
               sdata[tid] = sdata[tid+stride];
            }
         }
         __syncthreads();
      }

      // Store the output of the threads in this block
      if (threadIdx.x == 0) {
         max[bid] = sdata[0];
      }
   }
}

static __global__ void reductionMin(cudaReal* min, const cudaReal* in, int size)
{
   // number of blocks cut in two to avoid inactive initial threads
   int tid = threadIdx.x;
   int bid = blockIdx.x;
   int idx = bid * (blockDim.x*2) + tid;
   
   if (idx < size) {
      // Shared memory holding area
      extern __shared__ cudaReal sdata[];

      // Copy data into shared, wait for all threads to finish
      sdata[tid] = (in[idx] < in[idx+blockDim.x] ) ? in[idx] : in[idx + blockDim.x];
      
      __syncthreads();

      // Make comparisons across the block of data, each thread handling 
      // one comparison across two data points with strided indices before 
      // syncing with each other and then making further comparisons.
      for (int stride = blockDim.x/2; stride > 0; stride/=2) {
         if (tid < stride) {
            if (sdata[tid+stride] < sdata[tid]) {
               sdata[tid] = sdata[tid+stride];
            }
         }
         __syncthreads();
      }

      // Store the output of the threads in this block
      if (threadIdx.x == 0) {
         min[bid] = sdata[0];
      }
   }
}

static __global__ void reductionInnerProduct(cudaReal* innerprod, const cudaReal* a, const cudaReal* b int size)
{
   // number of blocks cut in two to avoid inactive initial threads
   int tid = threadIdx.x;
   int bid = blockIdx.x;
   int idx = bid * (blockDim.x*2) + tid;
   
   if (idx < size) {
      // Shared memory holding area
      extern __shared__ cudaReal sdata[];

      // Copy data into shared, wait for all threads to finish
      sdata[tid] = a[idx]*b[idx] + a[idx+blockDim.x]*b[idx+blockDim.x];
      
      __syncthreads();

      // Make comparisons across the block of data, each thread handling 
      // one comparison across two data points with strided indices before 
      // syncing with each other and then making further comparisons.
      for (int stride = blockDim.x/2; stride > 0; stride/=2) {
         if (tid < stride) {
            sdata[tid] += sdata[tid+stride]
         }
         __syncthreads();
      }

      // Store the output of the threads in this block
      if (threadIdx.x == 0) {
         innerprod[bid] = sdata[0];
      }
   }
}

template<unsigned int blockSize>
__global__ void deviceInnerProduct(cudaReal* c, const cudaReal* a,
   const cudaReal* b, int size) 
{
   //int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;

   //do all pointwise multiplication
   volatile extern __shared__ cudaReal cache[];
   cache[threadIdx.x] = a[startID] * b[startID];

   __syncthreads();

   if(blockSize >= 512) {
      if (threadIdx.x < 256){
         cache[threadIdx.x] += cache[threadIdx.x + 256];
      }
      __syncthreads();
   }
   if(blockSize >= 256) {
      if (threadIdx.x < 128){
         cache[threadIdx.x] += cache[threadIdx.x + 128];
      }
      __syncthreads();
   }
   if(blockSize >= 128) {
      if (threadIdx.x < 64){
         cache[threadIdx.x] += cache[threadIdx.x + 64];
      }
      __syncthreads();
   }

   if (threadIdx.x < 32) {
      if(blockSize >= 64) cache[threadIdx.x] += cache[threadIdx.x + 32];
      if(blockSize >= 32) cache[threadIdx.x] += cache[threadIdx.x + 16];
      if(blockSize >= 16) cache[threadIdx.x] += cache[threadIdx.x + 8];
      if(blockSize >= 8) cache[threadIdx.x] += cache[threadIdx.x + 4];
      if(blockSize >= 4) cache[threadIdx.x] += cache[threadIdx.x + 2];
      if(blockSize >= 2) cache[threadIdx.x] += cache[threadIdx.x + 1];
   }

   if (threadIdx.x == 0) {
      c[blockIdx.x] = cache[0];
   }
}

// --- Older Kernels --- //

// static __global__ void reductionOld(cudaReal* c, const cudaReal* a, int size) 
// {
//    // number of blocks used cut in two
//    int tid = threadIdx.x;
//    int bid = blockIdx.x;
//    int idx = bid * (blockDim.x*2) + tid;
   
   
//    volatile extern __shared__ cudaReal sdata[];
//    // global mem load combined with first operation
//    sdata[tid] = a[idx] + a[idx + blockDim.x];
   
//    __syncthreads();
   
//    // reduction
//    // data and block dimensions need to be a power of two
//    for (int stride = blockDim.x / 2; stride > 32; stride /= 2) {
//       if (tid < stride) {
//          sdata[tid] += sdata[tid + stride];
//       }
//       __syncthreads();
//    }
   
//    // unrolled. In final warp, synchronization is inherent 
//    if (tid < 32) {
      
//       sdata[tid] += sdata[tid + 32];
//       sdata[tid] += sdata[tid + 16];
//       sdata[tid] += sdata[tid + 8];
//       sdata[tid] += sdata[tid + 4];
//       sdata[tid] += sdata[tid + 2];
//       sdata[tid] += sdata[tid + 1];
      
//    }
   
//    if (tid == 0) {
//       c[bid] = sdata[0];
//    }
// }

#endif
