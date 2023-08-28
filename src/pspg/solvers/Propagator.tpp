#ifndef PSPG_PROPAGATOR_TPP
#define PSPG_PROPAGATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.h"
#include "Block.h"
#include <thrust/reduce.h>
#include "device_launch_parameters.h"
#include <cuda.h>
//#include <device_functions.h>
#include <thrust/count.h>
#include <pspg/math/GpuResources.h>
#include <pscf/mesh/Mesh.h>
//#include <Windows.h>

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
      if (isAllocated_) {
	      delete[] temp_;
	      cudaFree(d_temp_);
         cudaFree(qFields_d);
      }
   }  

   template <int D>
   void Propagator<D>::allocate(int ns, const Mesh<D>& mesh)
   {
      ns_ = ns;
      meshPtr_ = &mesh;

      ThreadGrid::setThreadsLogical(mesh.size());

      gpuErrchk(cudaMalloc((void**)&qFields_d, sizeof(cudaReal)* mesh.size() *
                 ns));
	   gpuErrchk(cudaMalloc((void**)&d_temp_, ThreadGrid::nBlocks() * sizeof(cudaReal)));
	   temp_ = new cudaReal[ThreadGrid::nBlocks()];
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

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(nx, nBlocks, nThreads);

      //qh[ix] = 1.0;
      //qFields_d points to the first element in gpu memory
      assignUniformReal<<<nBlocks, nThreads>>>(qFields_d, 1.0, nx);
      

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
         inPlacePointwiseMul<<<nBlocks, nThreads>>>(qFields_d, qt, nx);
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

      // GPU resources
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(nx, nBlocks, nThreads);

      // Initialize initial (head) field
      cudaReal* qh = qFields_d;
      // qh[i] = head[i];
      assignReal<<<nBlocks, nThreads>>>(qh, head, nx);

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
      
      Q = gpuInnerProduct(qh, qt, nx);
      Q /= double(nx);
      return Q;
   }

}
}
#endif
