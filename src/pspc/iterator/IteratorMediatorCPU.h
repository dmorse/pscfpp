#ifndef PSPC_ITERATOR_MEDIATOR_CPU_H
#define PSPC_ITERATOR_MEDIATOR_CPU_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <util/containers/DArray.h>
#include <pscf/iterator/IteratorMediator.h>

namespace Pscf {
namespace Pspc{

   using namespace Util;

   typedef DArray<double> FieldCPU;

   class IteratorMediatorCPU : public IteratorMediator<FieldCPU>
   {
   public:

      /// Constructor
      IteratorMediatorCPU(AbstractSystem& sys, Iterator<FieldCPU>& iter);

      /// Destructor
      ~IteratorMediatorCPU(); 

      /// Checks if the system has an initial guess
      bool hasInitialGuess();
      
      /// Calculates and returns the number of elements in the
      /// array to be iterated 
      int nElements();

      /// Gets a reference to the current state of the system
      void getCurrent(FieldCPU& curr);

      /// Runs calculation to evaluate function for fixed point.
      void evaluate();

      /// Gets residual values from system
      void getResidual(FieldCPU& resid);

      /// Updates the system with a passed in state of the iterator.
      void update(FieldCPU& newGuess);

   private:

      double scaleStress_ = 10.0;

   };

}
}
#endif