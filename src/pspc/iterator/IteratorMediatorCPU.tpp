#ifndef PSPC_ITERATOR_MEDIATOR_CPU_CPP
#define PSPC_ITERATOR_MEDIATOR_CPU_CPP

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

namespace Pscf {
namespace Pspc{

   using namespace Util;

   /// Constructor
   IteratorMediatorCPU::IteratorMediatorCPU(AbstractSystem& sys, Iterator<FieldCPU>& iter)
    : IteratorMediator<FieldCPU>(sys, iter)
   {}

   /// Destructor
   IteratorMediatorCPU::~IteratorMediatorCPU()
   {} 

   /// Checks if the system has an initial guess
   bool IteratorMediatorCPU::hasInitialGuess()
   {
      return system().hasWFields();
   }
   
   /// Calculates and returns the number of elements in the
   /// array to be iterated 
   int IteratorMediatorCPU::nElements()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().basis().nBasis();

      int nEle = nMonomer*nBasis;

      if (system().domain().isFlexible()) {
         nEle += system().unitCell().nParameter();
      }

      return nEle;
   }

   /// Gets a reference to the current state of the system
   void IteratorMediatorCPU::getCurrent(FieldCPU& curr)
   {
      // Straighten out fields into linear arrays

      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().basis().nBasis();
      const DArray<DArray<double>> * currSys = &system().wFields(); 
      
      
      for (int i = 0; i < nMonomer; i++) {
         for (int k = 0; k < nBasis; k++)
         {
            curr[i*nBasis+k] = (*currSys)[i][k];
         }
      }

      // If flexible unit cell, append them to the end of the current state
      if (system().domain().isFlexible()) {
         const int nParam = system().unitCell().nParameter();
         const FSArray<double,6> currParam = system().unitCell().parameters();

         for (int i = 0; i < nParam; i++) {
            curr[nMonomer*nBasis + i] = scaleStress_*currParam[i];
         }
      }

      return;
   }

   /// Runs calculation to evaluate function for fixed point.
   void IteratorMediatorCPU::evaluate()
   {
      system().compute();
      if (system().domain().isFlexible()) {
         system().mixture().computeStress();
      }
   }

   /// Gets residual values from system
   void IteratorMediatorCPU::getResidual(FieldCPU& resid)
   {
      const int n = nElements();
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().basis().nBasis();

      // Initialize residuals 
      for (int i = 0 ; i < n; ++i) {
         resid[i] = 0.0;
      }

      // Compute SCF residual vector elements
      for (int i = 0; i < nMonomer; ++i) {
         for (int j = 0; j < nMonomer; ++j) {
            for (int k = shift_; k < nBasis; ++k) {
               int idx = i*nBasis + k;
               resid[idx] +=
                  system().interaction().chi(i,j)*system().cField(j)[k] -
                  system().interaction().idemp(i,j)*system().wField(j)[k];
            }
         }
      }

      // If not canonical, account for incompressibility 
      if (!system().mixture().isCanonical()) {
         for (int i = 0; i < nMonomer; ++i) {
            resid[i*nBasis] -= 1.0/system().interaction().sum_inv();
         }
      } else {
         // otherwise explicitly set the residual value for the homogeneous components
         for (int i = 0; i < nMonomer; ++i) {
            resid[i*nBasis] = 0.0;
         }
      }

      // If variable unit cell, compute stress residuals
      if (system().domain().isFlexible()) {
         const int nParam = system().unitCell().nParameter();
         
         for (int i = 0; i < nParam ; i++) {
            resid[nMonomer*nBasis + i] = scaleStress_ * -1 * system().mixture().stress(i);
         }
      }

   }

   /// Updates the system with a passed in state of the iterator.
   void IteratorMediatorCPU::update(FieldCPU& newGuess)
   {
      // Convert back to field format
      const int nMonomer = system().mixture().nMonomer();
      const int nBasis = system().basis().nBasis();
      
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
      if (system().mixture().isCanonical()) {
         for (int i = 0; i < nMonomer; ++i) {
            wField[i][0] = 0.0; // initialize to 0
            for (int j = 0; j < nMonomer; ++j) {
               wField[i][0] += 
                  system().interaction().chi(i,j) * system().cField(j)[0];
            }
         }
      }
      system().setWBasis(wField);

      if (system().domain().isFlexible()) {
         FSArray<double, 6> parameters;
         const int nParam = system().unitCell().nParameter();

         for (int i = 0; i < nParam; i++) {
            parameters[i] = 1/scaleStress_ * newGuess[nMonomer*nBasis + i];
         }

         system().setUnitCell(parameters);
      }
      
   }

   /// Finalize????

}
}
#endif