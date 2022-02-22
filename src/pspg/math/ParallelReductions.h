#ifndef PARALLEL_REDUCTIONS_H
#define PARALLEL_REDUCTIONS_H

static __global__ void reductionSum(cudaReal* sum, const cudaReal* in, int size)
{
   // number of blocks cut in two to avoid inactive initial threads
   int tid = threadIdx.x;
   int bid = blockIdx.x;
   int idx = bid * (blockDim.x*2) + tid;
   
   if (idx+blockDim.x < size) {
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
      if (threadIdx.x == 0)
         sum[bid] = sdata[0];
   }
}

static __global__ void reductionMaxAbs(cudaReal* max, const cudaReal* in, int size) 
{
   // number of blocks cut in two to avoid inactive initial threads
   int tid = threadIdx.x;
   int bid = blockIdx.x;
   int idx = bid * (blockDim.x*2) + tid;
   
   if (idx+blockDim.x < size) {
      extern __shared__ cudaReal sdata[];
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
}

static __global__ void reductionMax(cudaReal* max, const cudaReal* in, int size)
{
   // number of blocks cut in two to avoid inactive initial threads
   int tid = threadIdx.x;
   int bid = blockIdx.x;
   int idx = bid * (blockDim.x*2) + tid;
   
   if (idx+blockDim.x < size) {
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
      if (threadIdx.x == 0)
         max[bid] = sdata[0];
   }
}

static __global__ void reductionMin(cudaReal* min, const cudaReal* in, int size)
{
   // number of blocks cut in two to avoid inactive initial threads
   int tid = threadIdx.x;
   int bid = blockIdx.x;
   int idx = bid * (blockDim.x*2) + tid;
   
   if (idx+blockDim.x < size) {
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
   
   if (idx+blockDim.x < size) {
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

#endif
