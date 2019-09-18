#ifndef GPU_RESOURCES_H
#define GPU_RESOURCES_H

extern int THREADS_PER_BLOCK;
extern int NUMBER_OF_BLOCKS;
#define SINGLE_PRECISION

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__);}
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
  if (code != cudaSuccess)
    {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if(abort) exit(code);
    }
}

static __global__ void helmholtzHelper(cufftReal* result, const cufftReal* composition,
   const cufftReal* pressure, float chi, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] = composition[i] * composition[i] / chi - pressure[i];
   }
}

static __global__ void reformField(cufftReal* Wa, cufftReal* Wb,
   const cufftReal* pressureF,const cufftReal* compositionF, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      Wa[i] = pressureF[i] + compositionF[i];
      Wb[i] = pressureF[i] - compositionF[i];
   }
}

static __global__ void mcftsStepHelper(cufftReal* result, const cufftReal* A2, const cufftReal* sqrtSq, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
          result[i] += A2[i] * sqrtSq[i];
   }
}
static __global__ void mcftsScale(cufftReal* result, cufftReal scale, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
          result[i] = result[i] * 2 * scale - scale;
   }
}

static __global__ void pointWiseAdd(cufftReal* result, const cufftReal* rhs, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] += rhs[i];
   }
}

static __global__ void subtractUniform(cufftReal* result, cufftReal rhs, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] -= rhs;
   }
}

static __global__ void pointWiseSubtract(cufftReal* result, const cufftReal* rhs, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] -= rhs[i];
   }
}

static __global__ void pointWiseSubtractFloat(cufftReal* result, const float rhs, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] -= rhs;
   }   
}

static __global__ void pointWiseBinarySubtract(const cufftReal* a, const cufftReal* b, cufftReal* result, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] = a[i] - b[i];
   }
}

static __global__ void pointWiseBinaryAdd(const cufftReal* a, const cufftReal* b, cufftReal* result, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] = a[i] + b[i];
   }
}

static __global__ void pointWiseAddScale(cufftReal* result, const cufftReal* rhs, float scale, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      result[i] += scale * rhs[i];
   }
}

//the 1 is a placeholder for dr
static __global__ void AmIsConvergedHelper(cufftReal* out, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   cufftReal temp;
   for (int i = startID; i < size; i += nThreads) {
      temp = (out[i] - 1) * (out[i] - 1) * 1;
      out[i] = temp;
   }
}

static __global__ void AmHelper(cufftReal* out, cufftReal* present, cufftReal* iPast, cufftReal* jPast, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      out[i] += (present[i] - iPast[i]) * (present[i] - jPast[i]);
   }
}

static __global__ void AmHelperVm(cufftReal* out, cufftReal* present, cufftReal* iPast, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for (int i = startID; i < size; i += nThreads) {
      out[i] += (present[i] - iPast[i]) * (present[i]);
   }
}

static __global__ void reduction(cufftReal* c, const cufftReal* a, int size) {
   //int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * (blockDim.x*2) + threadIdx.x;
   
   
   volatile extern __shared__ cufftReal cache[];
   //cufftReal temp = 0;
   //no need for loop here will be wrong.
   //for (int i = startID; i < size; i += nThreads) {
   //temp += a[startID];
   //}
   cache[threadIdx.x] = a[startID] + a[startID + blockDim.x];
   
   __syncthreads();
   
   //reduce operation
   //256/2 -- needs to be power of two
   for (int j = blockDim.x / 2; j > 32; j /= 2) {
      if (threadIdx.x < j) {
         cache[threadIdx.x] += cache[threadIdx.x + j];
      }
      __syncthreads();
   }
   
   if (threadIdx.x < 32) {
      
      cache[threadIdx.x] += cache[threadIdx.x + 32];
      cache[threadIdx.x] += cache[threadIdx.x + 16];
      cache[threadIdx.x] += cache[threadIdx.x + 8];
      cache[threadIdx.x] += cache[threadIdx.x + 4];
      cache[threadIdx.x] += cache[threadIdx.x + 2];
      cache[threadIdx.x] += cache[threadIdx.x + 1];
      
   }
   
   if (threadIdx.x == 0) {
      c[blockIdx.x] = cache[0];
   }
}

#endif
