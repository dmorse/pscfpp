#ifndef RP_ITERATOR_TPP
#define RP_ITERATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Iterator.h"
#include <prdc/environment/Environment.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   // Public functions

   /*
   * Default constructor.
   */
   template <int D, class T>
   Iterator<D,T>::Iterator()
    : isSymmetric_(false),
      isFlexible_(false),
      sysPtr_(nullptr)
   {  setClassName("Iterator"); }

   /*
   * Constructor.
   */
   template <int D, class T>
   Iterator<D,T>::Iterator(typename T::System& system)
    : isSymmetric_(false),
      isFlexible_(false),
      sysPtr_(&system)
   {  setClassName("Iterator"); }

   /*
   * Get the number of flexible lattice parameters.
   */
   template <int D, class T>
   int Iterator<D,T>::nFlexibleParams() const
   {
      UTIL_CHECK(sysPtr_);
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
   template <int D, class T>
   void Iterator<D,T>::setFlexibleParams(FSArray<bool,6> const & flexParams)
   {  
      UTIL_CHECK(sysPtr_);
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
   template <int D, class T>
   double Iterator<D,T>::stress(int paramId) const
   {
      // Parameter must be flexible to access the stress
      UTIL_CHECK(sysPtr_);
      UTIL_CHECK(flexibleParams_[paramId]);

      if (system().hasEnvironment()) {
         return system().environment().stress(paramId);
      } else {
         return system().mixture().stress(paramId);
      }
   }

   // Protected function

   /*
   * Set system.
   */
   template <int D, class T>
   void Iterator<D,T>::setSystem(typename T::System& system)
   {
      UTIL_CHECK(!sysPtr_);  
      sysPtr_ = &system; 
   }

} // namespace Rp
} // namespace Pscf
#endif
