#ifndef PSPG_ITERATOR_MEDIATOR_CUDA_TPP
#define PSPG_ITERATOR_MEDIATOR_CUDA_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include "IteratorMediatorCUDA.h"
#include <pscf/inter/ChiInteraction.h>

namespace Pscf {
namespace Pspg{

   using namespace Util;

   /// Constructor
   template <int D>
   IteratorMediatorCUDA<D>::IteratorMediatorCUDA(System<D>& sys)
    : sys_(&sys),
      iter_(0),
      temp_(0)
   {}

   /// Destructor
   template <int D>
   IteratorMediatorCUDA<D>::~IteratorMediatorCUDA()
   {
      if (iter_)
         delete iter_;

      if (temp_) {
         delete[] temp_;
         cudaFree(d_temp_);
      }
   }

   /// Set iterator pointer
   template <int D>
   void IteratorMediatorCUDA<D>::setIterator(Iterator<FieldCUDA>& iter)
   {
      iter_ = &iter;
      return;
   }

   template <int D>
   void IteratorMediatorCUDA<D>::setup()
   { iter_->setup(); }

   template <int D>
   int IteratorMediatorCUDA<D>::solve()
   { return iter_->solve(); }

   template <int D>
   bool IteratorMediatorCUDA<D>::hasInitialGuess()
   { return sys_->hasWFields(); }
   
   template <int D>
   int IteratorMediatorCUDA<D>::nElements()
   {
      const int nMonomer = sys_->mixture().nMonomer();
      const int nMesh = sys_->mesh().size();

      int nEle = nMonomer*nMesh;

      if (sys_->domain().isFlexible()) {
         nEle += sys_->unitCell().nParameter();
      }

      return nEle;
   }

   template <int D>
   void IteratorMediatorCUDA<D>::getCurrent(FieldCUDA& curr)
   {
      const int nMonomer = sys_->mixture().nMonomer();
      const int nMesh = sys_->mesh().size();
      const int n = nElements();

      // pointer to fields on system
      DArray<RDField<D>> * const currSys = &sys_->wFieldsRGrid();

      // Apply averaging convention from original GPU code -- set spatial average of
      // omega fields to 0 in the actual system.
      cudaReal average = 0;
      for (int i = 0; i < nMonomer; i++) {
         // Find current spatial average
         average += findAverage((*currSys)[i].cDField(), nMesh);
      } 
      average /= nMonomer;
      // Subtract average from each field, setting average to zero
      for (int i = 0; i < nMonomer; i++) {
         subtractUniform<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>((*currSys)[i].cDField(), average, nMesh);
      } 
      
      // Initialize current field values values to zeroes.
      gpuErrchk( cudaMemset(curr.cDField(), 0, n*sizeof(cudaReal)) );

      // loop to unfold the system fields and store them in one long array
      for (int i = 0; i < nMonomer; i++) {
         assignReal<<<NUMBER_OF_BLOCKS,THREADS_PER_BLOCK>>>(curr.cDField() + i*nMesh, (*currSys)[i].cDField(), nMesh);
      }

      // if flexible unit cell, also store unit cell parameters
      if (sys_->domain().isFlexible()) {
         const int nParam = sys_->unitCell().nParameter();
         const double scaleStress = sys_->domain().scaleStress();
         const FSArray<double,6> currParam = sys_->unitCell().parameters();
         // convert into a cudaReal array
         cudaReal* temp = new cudaReal[nParam];
         for (int k = 0; k < nParam; k++) 
               temp[k] = (cudaReal)scaleStress*currParam[k];
         
         // copy paramters to the end of the curr array
         gpuErrchk(cudaMemcpy(curr.cDField() + nMonomer*nMesh, temp, nParam*sizeof(cudaReal), cudaMemcpyHostToDevice));
         delete[] temp;
      }
   }

   template <int D>
   void IteratorMediatorCUDA<D>::evaluate()
   {
      // Solve MDEs for current omega field
      sys_->compute();
      // Compute stress if done
      if (sys_->domain().isFlexible()) {
         sys_->mixture().computeStress(sys_->wavelist());
      }
   }

   template <int D>
   void IteratorMediatorCUDA<D>::getResidual(FieldCUDA& resid)
   {
      const int n = nElements();
      const int nMonomer = sys_->mixture().nMonomer();
      const int nMesh = sys_->mesh().size();
      
      // Initialize residuals to zero
      assignUniformReal<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(resid.cDField(), 0, n);

      // Compute SCF residuals
      for (int i = 0; i < nMonomer; i++) {
         int startIdx = i*nMesh;
         for (int j = 0; j < nMonomer; j++) {
            pointWiseAddScale<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(resid.cDField() + startIdx,
                                                                            sys_->cFieldRGrid(j).cDField(),
                                                                            sys_->interaction().chi(i, j),
                                                                            nMesh);
            pointWiseAddScale<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(resid.cDField() + startIdx,
                                                                            sys_->wFieldRGrid(j).cDField(),
                                                                            -sys_->interaction().idemp(i, j),
                                                                            nMesh);
         }
      }
      

      // Always do this subtraction, regardless of whether canonical. From original GPU code.
      cudaReal factor = 1/(cudaReal)sys_->interaction().sum_inv();
      for (int i = 0; i < nMonomer; ++i) {
         
         subtractUniform<<<NUMBER_OF_BLOCKS, THREADS_PER_BLOCK>>>(resid.cDField() + i*nMesh,
                                                                           factor, 
                                                                           nMesh);
      }

      // If variable unit cell, compute stress residuals
      if (sys_->domain().isFlexible()) {
         const int nParam = sys_->unitCell().nParameter();
         const double scaleStress = sys_->domain().scaleStress();
         cudaReal* stress = new cudaReal[nParam];

         for (int i = 0; i < nParam; i++) {
            stress[i] = (cudaReal)(-1*scaleStress*sys_->mixture().stress(i));
         }

         gpuErrchk(cudaMemcpy(resid.cDField()+nMonomer*nMesh, stress, nParam*sizeof(cudaReal), cudaMemcpyHostToDevice));
      }
   }

   template <int D>
   void IteratorMediatorCUDA<D>::update(FieldCUDA& newGuess)
   {
      const int nMonomer = sys_->mixture().nMonomer();
      const int nMesh = sys_->mesh().size();

      // pointer to fields on system
      DArray<RDField<D>> * sysfields = &sys_->wFieldsRGrid();

      // copy over grid points
      for (int i = 0; i < nMonomer; i++) {
         assignReal<<<NUMBER_OF_BLOCKS,THREADS_PER_BLOCK>>>
               ((*sysfields)[i].cDField(), newGuess.cDField() + i*nMesh, nMesh);
      }

      // if flexible unit cell, update parameters well
      if (sys_->domain().isFlexible()) {
         FSArray<double,6> parameters;
         const int nParam = sys_->unitCell().nParameter();
         const double scaleStress = sys_->domain().scaleStress();
         cudaReal* temp = new cudaReal[nParam];

         gpuErrchk(cudaMemcpy(temp, newGuess.cDField() + nMonomer*nMesh, nParam*sizeof(cudaReal), cudaMemcpyDeviceToHost));
         for (int i = 0; i < nParam; i++) {
            parameters.append(1/scaleStress * (double)temp[i]);
         }

         sys_->unitCell().setParameters(parameters);            
         sys_->mixture().setupUnitCell(sys_->unitCell(), sys_->wavelist());
         sys_->wavelist().computedKSq(sys_->unitCell());

         delete[] temp;
      }

   }

   // --- Private member functions that are specific to this implementation --- 

   template<int D> 
   cudaReal IteratorMediatorCUDA<D>::findAverage(cudaReal * const field, int n) 
   {
      cudaReal average = gpuSum(field, n)/n;

      return average;
   }

   
   // Finalize???? Could have the iterator call this, and have it be by default empty
   // as a place to put things that may need to be done at the end of an iteration. Or could have the system
   // take care of this since it doesn't have much to do with the iterator.
}
}
#endif