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
      iter_(0)
   {}

   /// Destructor
   template <int D>
   IteratorMediatorCUDA<D>::~IteratorMediatorCUDA()
   {
      if (iter_)
         delete iter_;
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
      int nEle = nData();


      // If dataset is not a power of two, decide what to do. We can either have the GPU
      // handle the excess with a bunch of zeroes to make it a power of two, or we can have
      // the CPU handle the excess by manually processing the excess.
      if ( (n & (n - 1)) != 0 ) {
         double pow = log2(nEle);
         // distance to nearest higher power of two
         int nExCeil = 2**ceil(pow) - nEle;
         // distance to nearest lower power of two
         int nExFloor = nEle - 2**floor(pow);
         // A tunable parameter to indicate the relative
         // cost of a GPU data point versus a CPU data point.
         // I expect this to be less than 1  
         int weightGPU = 1;
         // if additional data points exceed the weighted 
         // number of additional points to be added to make it 
         // a power of two, then force it to a power of two
         if (nExFloor > weightGPU*nExCeil)
            nEle = 2**ceil(pow);
         // otherwise, leave nEle the same. Excess elements will be handled by the CPU, as 
         // deined in AmStrategyCUDA
      }
        
      return nEle;
   }

   template <int D>
   void IteratorMediatorCUDA<D>::getCurrent(FieldCUDA& curr)
   {
      const int nMonomer = sys_->mixture().nMonomer();
      const int nMesh = sys_->mesh().size();

      // Initialize values to zeroes. This is important especially in the case of 
      // the number of elements being larger to make the system a power of
      cudaMemset(curr.cDField(), 0, nElements()*sizeof(cudaReal));

      // pointer to fields on system
      const DArray<RDField<D>> * currSys = &sys_->wFieldsRGrid();

      // loop to unfold the system fields and store them in one long array
      for (int i = 0; i < nMonomer; i++) {
         assignReal<<<NUMBER_OF_BLOCKS,THREADS_PER_BLOCK>>>(curr + i*nMesh, (*currSys)[i].cDField(), nMesh);
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
         
         // copy paramters the end of the curr array
         cudaMemcpy(curr.cDField() + nMonomer*nMesh, temp, nParam*sizeof(cudaReal), cudaMemcpyHostToDevice)
         delete temp;
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
      // calculate residuals. why is this done here and not in strategy? bc we need to compare fields at the system level
   }

   template <int D>
   void IteratorMediatorCUDA<D>::update(FieldCUDA& newGuess)
   {
      // push new fields to system.

      const int nMonomer = sys_->mixture().nMonomer();
      const int nMesh = sys_->mesh().size();

      // pointer to fields on system
      DArray<RDField<D>> * sysfields = &sys_->wFieldsRGrid();
      
      // copy over grid points
      for (int i = 0; i < nMonomer; i++) {
         assignReal<<<NUMBER_OF_BLOCKS,THREADS_PER_BLOCK>>>
               ((*sysfields)[i], newGuess.cDField() + i*nMesh, nMesh);
      }

      // if flexible unit cell, update parameters well
      if (sys_->domain().isFlexible()) {
         FSArray<double,6> newParam;
         cudaReal* temp = new cudaReal[nParam];

         cudaMemcpy(temp, newGuess.cDField() + nMonomer*nMesh, nParam*sizeof(cudaReal), cudaMemcpyDeviceToHost);
         for (int i = 0; i < nParam; i++) {
            newParam.append((double)temp[i]);
         }

         systemPtr_->unitCell().setParameters(newParam);            
         systemPtr_->mixture().setupUnitCell(systemPtr_->unitCell(), systemPtr_->wavelist());
         systemPtr_->wavelist().computedKSq(systemPtr_->unitCell());

         delete temp;
      }

   }
   
   template <int D>
   int IteratorMediatorCUDA<D>::nData() 
   {
      const int nMonomer = sys_->mixture().nMonomer();
      const int nMesh = sys_->mesh().size();

      int nEle = nMonomer*nMesh;

      if (sys_->domain().isFlexible()) {
         nEle += sys_->unitCell().nParameter();
      }

      return nEle;
   }

   // Finalize???? Could have the iterator call this, and have it be by default empty
   // as a place to put things that may need to be done at the end of an iteration. Or could have the system
   // take care of this since it doesn't have much to do with the iterator.
}
}
#endif