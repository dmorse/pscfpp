#ifndef PSPG_PROPAGATOR_TPP
#define PSPG_PROPAGATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.h"
#include "Block.h"
#include <thrust/reduce.h>
#include "device_launch_parameters.h"
#include <cuda.h>
//#include <device_functions.h>
#include <thrust/count.h>
#include <pspg/GpuResources.h>
#include <pscf/mesh/Mesh.h>
//#include <Windows.h>

__global__ 
void assignUniformReal(cudaReal* result, cudaReal uniform, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < size; i += nThreads) {
      result[i] = uniform;
   }
}

__global__ 
void assignReal(cudaReal* result, const cudaReal* rhs, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < size; i += nThreads) {
      result[i] = rhs[i];
   }
}

__global__ 
void inPlacePointwiseMul(cudaReal* a, const cudaReal* b, int size) {
   int nThreads = blockDim.x * gridDim.x;
   int startID = blockIdx.x * blockDim.x + threadIdx.x;
   for(int i = startID; i < size; i += nThreads) {
      a[i] *= b[i];
   }
}

#if 0
template<unsigned int blockSize>
__global__ void deviceInnerProduct(cudaReal* c, const cudaReal* a,
	const cudaReal* b, int size) {
        //int nThreads = blockDim.x * gridDim.x;
	int startID = blockIdx.x * blockDim.x + threadIdx.x;

	//do all pointwise multiplication
	volatile extern __shared__ cudaReal cache[];
	cudaReal temp = 0;
	//no need for loop here will be wrong.
	//for (int i = startID; i < size; i += nThreads) {
	temp += a[startID] * b[startID];
	//}
	cache[threadIdx.x] = temp;

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
	//reduce operation
	//256/2 -- needs to be power of two
	//for (int j = blockDim.x / 2; j > 32; j /= 2) {
	//	if (threadIdx.x < j) {
	//		cache[threadIdx.x] += cache[threadIdx.x + j];
	//	}
	//	__syncthreads();
	//}
	

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

namespace Pscf {
namespace Pspg {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   Propagator<D>::Propagator()
    : blockPtr_(0),
      meshPtr_(0),
      ns_(0),
      temp_(0),
      isAllocated_(false)
   {
   }

   /*
   * Destructor.
   */
   template <int D>
   Propagator<D>::~Propagator()
   {
      if (temp_) {
	      delete[] temp_;
	      cudaFree(d_temp_);
      }
   }  

   template <int D>
   void Propagator<D>::allocate(int ns, const Mesh<D>& mesh)
   {
      ns_ = ns;
      meshPtr_ = &mesh;

      cudaMalloc((void**)&qFields_d, sizeof(cudaReal)* mesh.size() *
                 ns);
	   cudaMalloc((void**)&d_temp_, NUMBER_OF_BLOCKS * sizeof(cudaReal));
	   temp_ = new cudaReal[NUMBER_OF_BLOCKS];
      isAllocated_ = true;
   }

   /*
   * Compute initial head QField from final tail QFields of sources.
   */
   template <int D>
   void Propagator<D>::computeHead()
   {

      // Reference to head of this propagator
      //QField& qh = qFields_[0];

      // Initialize qh field to 1.0 at all grid points
      int nx = meshPtr_->size();

      //qh[ix] = 1.0;
      //qFields_d points to the first float in gpu memory
      assignUniformReal<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(qFields_d, 1.0, nx);
      

      // Pointwise multiply tail QFields of all sources
      // this could be slow with many sources. Should launch 1 kernel for the whole
      // function of computeHead
      const cudaReal* qt;
      for (int is = 0; is < nSource(); ++is) {
         if (!source(is).isSolved()) {
            UTIL_THROW("Source not solved in computeHead");
         }
         //need to modify tail to give the total_size - mesh_size pointer
         qt = source(is).tail();

         //qh[ix] *= qt[ix];
         inPlacePointwiseMul<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>
            (qFields_d, qt, nx);
      }
   }

   /*
   * Solve the modified diffusion equation for this block.
   */
   template <int D>
   void Propagator<D>::solve()
   {
      UTIL_CHECK(isAllocated());
      computeHead();
      // Setup solver and solve
      block().setupFFT();
      //cudaReal* qf;
      //qf = new cudaReal;
      
      int currentIdx;
      for (int iStep = 0; iStep < ns_ - 1; ++iStep) {
         currentIdx = iStep * meshPtr_->size();
         //block has to learn to deal with the cudaReal
         block().step(qFields_d + currentIdx, qFields_d + currentIdx + meshPtr_->size());
      }
	  //delete qf;
      setIsSolved(true);
   }

   /*
   * Solve the modified diffusion equation with specified initial field.
   */
   template <int D>
   void Propagator<D>::solve(const cudaReal * head)
   {
      int nx = meshPtr_->size();

      // Initialize initial (head) field
      cudaReal* qh = qFields_d;
      // qh[i] = head[i];
      assignReal<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(qh, head, nx);

      // Setup solver and solve
      int currentIdx;
      for (int iStep = 0; iStep < ns_ - 1; ++iStep) {
         currentIdx = iStep * nx;
         block().step(qFields_d + currentIdx, qFields_d + currentIdx + nx);
      }
      setIsSolved(true);
   }

   /*
   * Integrate to calculate monomer concentration for this block
   */
   template <int D>
   double Propagator<D>::computeQ()
   {
      // Preconditions
      if (!isSolved()) {
         UTIL_THROW("Propagator is not solved.");
      }
      if (!hasPartner()) {
         UTIL_THROW("Propagator has no partner set.");
      }
      if (!partner().isSolved()) {
         UTIL_THROW("Partner propagator is not solved");
      }
      const cudaReal * qh = head();
      const cudaReal * qt = partner().tail();
      int nx = meshPtr_->size();

      // Take inner product of head and partner tail fields
      // cannot reduce assuming one propagator, qh == 1
      // polymers are divided into blocks midway through
      double Q = 0; 
      
      Q = innerProduct(qh, qt, nx);
      Q /= double(nx);
      return Q;
   }

   template <int D>
   cudaReal Propagator<D>::innerProduct(const cudaReal* a, const cudaReal* b, int size) {
	   
     switch(THREADS_PER_BLOCK){
     case 512:
       deviceInnerProduct<512><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cudaReal)>>>(d_temp_, a, b, size);
       break;
     case 256:
       deviceInnerProduct<256><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cudaReal)>>>(d_temp_, a, b, size);
       break;
     case 128:
       deviceInnerProduct<128><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cudaReal)>>>(d_temp_, a, b, size);
       break;
     case 64:
       deviceInnerProduct<64><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cudaReal)>>>(d_temp_, a, b, size);
       break;
     case 32:
       deviceInnerProduct<32><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cudaReal)>>>(d_temp_, a, b, size);
       break;
     case 16:
       deviceInnerProduct<16><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cudaReal)>>>(d_temp_, a, b, size);
       break;
     case 8:
       deviceInnerProduct<8><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cudaReal)>>>(d_temp_, a, b, size);
       break;
     case 4:
       deviceInnerProduct<4><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cudaReal)>>>(d_temp_, a, b, size);
       break;
     case 2:
       deviceInnerProduct<2><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cudaReal)>>>(d_temp_, a, b, size);
       break;
     case 1:
       deviceInnerProduct<1><<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK, THREADS_PER_BLOCK * sizeof(cudaReal)>>>(d_temp_, a, b, size);
       break;
     }
     
	   cudaMemcpy(temp_, d_temp_, NUMBER_OF_BLOCKS * sizeof(cudaReal), cudaMemcpyDeviceToHost);
	   cudaReal final = 0;
	   cudaReal c = 0;
	   //use kahan summation
	   for(int i = 0; i < NUMBER_OF_BLOCKS; ++i) {
		   cudaReal y = temp_[i] - c;
		   cudaReal t = final + y;
		   c = (t - final) - y;
		   final = t;
	   }
	   
	   return final;
   }

}
}
#endif
