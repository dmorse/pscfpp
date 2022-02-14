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

static __global__ void pointWiseAddScale(cudaReal* result, const cudaReal* rhs, float scale, int size)
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] += scale * rhs[i];
   }
}

//the 1 is a placeholder for dr
static __global__ void AmIsConvergedHelper(cudaReal* out, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   cudaReal temp;
   for (int i = startID; i < size; i += nThreads) {
      temp = (out[i] - 1) * (out[i] - 1) * 1;
      out[i] = temp;
   }
}

static __global__ void AmHelper(cudaReal* out, cudaReal* present, cudaReal* iPast, cudaReal* jPast, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      out[i] += (present[i] - iPast[i]) * (present[i] - jPast[i]);
   }
}

static __global__ void AmHelperVm(cudaReal* out, cudaReal* present, cudaReal* iPast, int size) 
{
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      out[i] += (present[i] - iPast[i]) * (present[i]);
   }
}

static __global__ void reduction(cudaReal* c, const cudaReal* a, int size) 
{
   // number of blocks used cut in two
   int tid = threadIdx.x;
   int bid = blockIdx.x;
   int idx = bid * (blockDim.x*2) + tid;
   
   
   volatile extern __shared__ cudaReal sdata[];
   // global mem load combined with first operation
   sdata[tid] = a[idx] + a[idx + blockDim.x];
   
   __syncthreads();
   
   // reduction
   // data and block dimensions need to be a power of two
   for (int stride = blockDim.x / 2; stride > 32; stride /= 2) {
      if (tid < stride) {
         sdata[tid] += sdata[tid + stride];
      }
      __syncthreads();
   }
   
   // unrolled. In final warp, synchronization is inherent 
   if (tid < 32) {
      
      sdata[tid] += sdata[tid + 32];
      sdata[tid] += sdata[tid + 16];
      sdata[tid] += sdata[tid + 8];
      sdata[tid] += sdata[tid + 4];
      sdata[tid] += sdata[tid + 2];
      sdata[tid] += sdata[tid + 1];
      
   }
   
   if (tid == 0) {
      c[bid] = sdata[0];
   }
}

template <typename T>
static __global__ void reduction(T* d_out, const T* d_in, T reduceFunc(T), int size) 
{
   // number of blocks cut in two to avoid inactive initial threads
   int tid = threadIdx.x;
   int bid = blockIdx.x;
   int idx = bid * (blockDim.x*2) + tid;
   
   
   volatile extern __shared__ T sdata[];
   // global mem load combined with first operation
   sdata[tid] = reduceFunc( d_in[idx], d_in[idx + blockDim.x] );
   
   __syncthreads();
   
   // reduction
   // data and block dimensions need to be a power of two
   for (int stride = blockDim.x / 2; stride > 32; stride /= 2) {
      if (tid < stride) {
         sdata[tid] = reduceFunc( sdata[tid], sdata[tid + stride] );
      }
      __syncthreads();
   }
   
   // unrolled. In final warp, synchronization is inherent.
   if (tid < 32) warpReduce(sdata, reduceFunc, tid);
   
   if (tid == 0) {
      c[bid] = sdata[0];
   }
}

template <typename T>
static __device__ void warpReduce(volatile T* sdata, T reduceFunc(T), int tid)
{
   sdata[tid] = reduceFunc( sdata[tid], sdata[tid + 32] );
   sdata[tid] = reduceFunc( sdata[tid], sdata[tid + 16] );
   sdata[tid] = reduceFunc( sdata[tid], sdata[tid +  8]  );
   sdata[tid] = reduceFunc( sdata[tid], sdata[tid +  4]  );
   sdata[tid] = reduceFunc( sdata[tid], sdata[tid +  2]  );
   sdata[tid] = reduceFunc( sdata[tid], sdata[tid +  1]  );
}

static __global__ void reductionMax(cudaReal* d_max, cudaReal* in, int size)
{
   int idx = blockIdx.x * blockDim.x + threadIdx.x;
   
   if (idx < size) {
      // Shared memory holding area
      extern __shared__ cudaReal sharedData[];

      // Copy data into shared, wait for all threads to finish
      sharedData[threadIdx.x] = in[idx];
      
      __syncthreads();

      // Make comparisons across the block of data, each thread handling 
      // one comparison across two data points with strided indices before 
      // syncing with each other and then making further comparisons.
      for (int stride = blockDim.x/2; stride > 0; stride/=2) {
         if (threadIdx.x < stride) {
            if (sharedData[threadIdx.x+stride] > sharedData[threadIdx.x]) {
               sharedData[threadIdx.x] = sharedData[threadIdx.x+stride];
            }
         }
         

         __syncthreads();
      }

      // Store the output of the threads in this block
      if (threadIdx.x == 0) {
         d_max[blockIdx.x] = sharedData[0];
      }
   }

}

static __global__ void reductionMin(cudaReal* d_min, cudaReal* in, int size)
{
   int idx = blockIdx.x * blockDim.x + threadIdx.x;
   
   if (idx < size) {
      // Shared memory holding area
      extern __shared__ cudaReal sharedData[];

      // Copy data into shared, wait for all threads to finish
      sharedData[threadIdx.x] = in[idx];
      
      __syncthreads();

      // Make comparisons across the block of data, each thread handling 
      // one comparison across two data points with strided indices before 
      // syncing with each other and then making further comparisons.
      for (int stride = blockDim.x/2; stride > 0; stride/=2) {
         if (threadIdx.x < stride) {
            if (sharedData[threadIdx.x+stride] < sharedData[threadIdx.x]) {
               sharedData[threadIdx.x] = sharedData[threadIdx.x+stride];
            }
         }
         

         __syncthreads();
      }

      // Store the output of the threads in this block
      if (threadIdx.x == 0) {
         d_max[blockIdx.x] = sharedData[0];
      }
   }

}

#endif
