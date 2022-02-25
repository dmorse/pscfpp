#ifndef PSCF_ITERATOR_MEDIATOR_H
#define PSCF_ITERATOR_MEDIATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <pscf/AbstractSystem.h>
#include <pscf/iterator/Iterator.h>

namespace Pscf {

   // Forward declarations
   class AbstractSystem;
   template <typename T> class Iterator;

   using namespace Util;

   template <typename T>
   class IteratorMediator 
   {
   public:

      /// Constructor
      IteratorMediator() {};

      /// Destructor
      virtual ~IteratorMediator() {}; 

      /// Set iterator pointer
      virtual void setIterator(Iterator<T>& iter) = 0;
      
      virtual void setup() = 0;

      virtual int solve() = 0;

      /// Checks if the system has an initial guess
      virtual bool hasInitialGuess() = 0;
      
      /// Calculates and returns the number of elements in the
      /// array to be iterated 
      virtual int nElements() = 0;

      /// Gets a reference to the current state of the system
      virtual void getCurrent(T& curr) = 0;

      /// Runs calculation to evaluate function for fixed point.
      virtual void evaluate() = 0;

      /// Gets residual values from system
      virtual void getResidual(T& resid) = 0;

      /// Updates the system with a passed in state of the iterator.
      virtual void update(T& newGuess) = 0;

   };

}
#endif