#ifndef PSPC_ITERATOR_MEDIATOR_CPU_TPP
#define PSPC_ITERATOR_MEDIATOR_CPU_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <util/containers/DArray.h>
#include <util/containers/FSArray.h>
#include "IteratorMediatorCPU.h"
#include <pscf/inter/ChiInteraction.h>

namespace Pscf {
namespace Pspc{

   using namespace Util;

   /// Constructor
   template <int D>
   IteratorMediatorCPU<D>::IteratorMediatorCPU(System<D>& sys, Pscf::Iterator<FieldCPU>& iter)
    : sys_(&sys),
      iter_(&iter),
      scaleStress_(10.0)
   {}

   /// Destructor
   template <int D>
   IteratorMediatorCPU<D>::~IteratorMediatorCPU()
   {}

   template <int D>
   void IteratorMediatorCPU<D>::setup()
   { iter_->setup(); }

   template <int D>
   int IteratorMediatorCPU<D>::solve()
   { return iter_->solve(); }

   template <int D>
   bool IteratorMediatorCPU<D>::hasInitialGuess()
   {
      return sys_->hasWFields();
   }
   
   template <int D>
   int IteratorMediatorCPU<D>::nElements()
   {
      const int nMonomer = sys_->mixture().nMonomer();
      const int nBasis = sys_->basis().nBasis();

      int nEle = nMonomer*nBasis;

      if (sys_->domain().isFlexible()) {
         nEle += sys_->unitCell().nParameter();
      }

      return nEle;
   }

   template <int D>
   void IteratorMediatorCPU<D>::getCurrent(FieldCPU& curr)
   {
      // Straighten out fields into linear arrays

      const int nMonomer = sys_->mixture().nMonomer();
      const int nBasis = sys_->basis().nBasis();
      const DArray<DArray<double>> * currSys = &sys_->wFields(); 
      
      
      for (int i = 0; i < nMonomer; i++) {
         for (int k = 0; k < nBasis; k++)
         {
            curr[i*nBasis+k] = (*currSys)[i][k];
         }
      }

      if (sys_->domain().isFlexible()) {
         const int nParam = sys_->unitCell().nParameter();
         const FSArray<double,6> currParam = sys_->unitCell().parameters();

         for (int i = 0; i < nParam; i++) {
            curr[nMonomer*nBasis + i] = scaleStress_*currParam[i];
         }
      }

      return;
   }

   template <int D>
   void IteratorMediatorCPU<D>::evaluate()
   {
      sys_->compute();
      if (sys_->domain().isFlexible()) {
         sys_->mixture().computeStress();
      }
   }

   template <int D>
   void IteratorMediatorCPU<D>::getResidual(FieldCPU& resid)
   {
      const int n = nElements();
      const int nMonomer = sys_->mixture().nMonomer();
      const int nBasis = sys_->basis().nBasis();

      // Initialize residuals 
      for (int i = 0 ; i < n; ++i) {
         resid[i] = 0.0;
      }

      // Compute SCF residual vector elements
      for (int i = 0; i < nMonomer; ++i) {
         for (int j = 0; j < nMonomer; ++j) {
            for (int k = 0; k < nBasis; ++k) {
               int idx = i*nBasis + k;
               resid[idx] +=
                  sys_->interaction().chi(i,j)*sys_->cField(j)[k] -
                  sys_->interaction().idemp(i,j)*sys_->wField(j)[k];
            }
         }
      }

      // If not canonical, account for incompressibility 
      if (!sys_->mixture().isCanonical()) {
         for (int i = 0; i < nMonomer; ++i) {
            resid[i*nBasis] -= 1.0/sys_->interaction().sum_inv();
         }
      } else {
         // otherwise explicitly set the residual value for the homogeneous components
         for (int i = 0; i < nMonomer; ++i) {
            resid[i*nBasis] = 0.0;
         }
      }

      // If variable unit cell, compute stress residuals
      if (sys_->domain().isFlexible()) {
         const int nParam = sys_->unitCell().nParameter();
         
         for (int i = 0; i < nParam ; i++) {
            resid[nMonomer*nBasis + i] = scaleStress_ * -1 * sys_->mixture().stress(i);
         }
      }

   }

   template <int D>
   void IteratorMediatorCPU<D>::update(FieldCPU& newGuess)
   {
      // Convert back to field format
      const int nMonomer = sys_->mixture().nMonomer();
      const int nBasis = sys_->basis().nBasis();
      
      DArray< DArray<double> > wField;
      wField.allocate(nMonomer);
      
      // Restructure in format of monomers, basis functions
      for (int i = 0; i < nMonomer; i++) {
         wField[i].allocate(nBasis);
         for (int k = 0; k < nBasis; k++)
         {
            wField[i][k] = newGuess[i*nBasis + k];
         }
      }
      // Manually and explicitly set homogeneous components of field if canonical
      if (sys_->mixture().isCanonical()) {
         for (int i = 0; i < nMonomer; ++i) {
            wField[i][0] = 0.0; // initialize to 0
            for (int j = 0; j < nMonomer; ++j) {
               wField[i][0] += 
                  sys_->interaction().chi(i,j) * sys_->cField(j)[0];
            }
         }
      }
      sys_->setWBasis(wField);

      if (sys_->domain().isFlexible()) {
         FSArray<double, 6> parameters;
         const int nParam = sys_->unitCell().nParameter();

         for (int i = 0; i < nParam; i++) {
            parameters[i] = 1/scaleStress_ * newGuess[nMonomer*nBasis + i];
         }

         sys_->setUnitCell(parameters);
      }
      
   }

   /// Finalize????

}
}
#endif