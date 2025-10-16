#ifndef RPG_ITERATOR_TPP
#define RPG_ITERATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"

#include <prdc/environment/Environment.h>
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>

namespace Pscf {
namespace Rpg {

   template <int D> class System;

   using namespace Util;
   using namespace Prdc;

   /*
   * Default constructor.
   */
   template<int D>
   Iterator<D>::Iterator()
    : isSymmetric_(false),
      isFlexible_(false),
      sysPtr_(nullptr)
   {  setClassName("Iterator"); }

   /*
   * Constructor.
   */
   template<int D>
   Iterator<D>::Iterator(System<D>& system)
    : isSymmetric_(false),
      isFlexible_(false),
      sysPtr_(&system)
   {  setClassName("Iterator"); }

   /*
   * Destructor.
   */
   template<int D>
   Iterator<D>::~Iterator()
   {}

   /*
   * Does this iterator use a symmetry-adapted Fourier basis?
   */
   template<int D>
   inline bool Iterator<D>::isSymmetric() const
   {  return isSymmetric_; }

   /*
   * Is the unit cell flexible (true) or rigid (false) ?
   */
   template<int D>
   inline bool Iterator<D>::isFlexible() const
   {  return isFlexible_; }

   /*
   * Get the array indicating which lattice parameters are flexible.
   */
   template<int D>
   inline FSArray<bool,6> Iterator<D>::flexibleParams() const
   {  return flexibleParams_; }

   /*
   * Get the number of flexible lattice parameters.
   */
   template <int D>
   int Iterator<D>::nFlexibleParams() const
   {
      UTIL_CHECK(flexibleParams_.size() == 
                                system().domain().unitCell().nParameter());
      int nFlexParams = 0;
      for (int i = 0; i < flexibleParams_.size(); i++) {
         if (flexibleParams_[i]) nFlexParams++;
      }
      return nFlexParams;
   }

   /*
   * Set the array indicating which lattice parameters are flexible.
   */
   template <int D>
   void Iterator<D>::setFlexibleParams(FSArray<bool,6> const & flexParams)
   {  
      flexibleParams_ = flexParams; 
      if (nFlexibleParams() == 0) {
         isFlexible_ = false;
      } else {
         isFlexible_ = true;
      }
   }

   /*
   * Return the stress used by this Iterator, for one lattice parameter.
   */
   template <int D>
   double Iterator<D>::stress(int paramId) const
   {
      // Parameter must be flexible to access the stress
      UTIL_CHECK(flexibleParams_[paramId]);

      if (system().hasEnvironment()) {
         return system().environment().stress(paramId);
      } else {
         return system().mixture().stress(paramId);
      }
   }

} // namespace Rpg
} // namespace Pscf
#endif
