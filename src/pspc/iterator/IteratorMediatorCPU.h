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
#include <pspc/System.h>

namespace Pscf {
namespace Pspc{

   using namespace Util;

   typedef DArray<double> FieldCPU;

   template <int D>
   class IteratorMediatorCPU : public IteratorMediator<FieldCPU>
   {
   public:

      /// Constructor
      IteratorMediatorCPU(System<D>& sys);

      /// Destructor
      ~IteratorMediatorCPU(); 

      /// Set iterator pointer
      void setIterator(Iterator<FieldCPU>& iter); 

      /// Set up the iterator.
      void setup();

      /// Instructor the iterator to solve the fixed point problem.
      int solve();

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

      /// Outputs relevant system details to the iteration log
      void outputToLog();

   private:

      // pointer to system
      System<D>* sys_;

      // pointer to iterator
      Pscf::Iterator<FieldCPU>* iter_; 

   };

   #ifndef PSPC_ITERATOR_MEDIATOR_CPU_TPP
   // Suppress implicit instantiation
   extern template class IteratorMediatorCPU<1>;
   extern template class IteratorMediatorCPU<2>;
   extern template class IteratorMediatorCPU<3>;
   #endif

}
}
#endif