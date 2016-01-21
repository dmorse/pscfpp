#ifndef FD1D_SYSTEM_H
#define FD1D_SYSTEM_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"
#include "Iterator.h"
#include <util/param/ParamComposite.h>

namespace Pscf {
namespace Fd1d
{

   using namespace Util;

   class System : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      System();

      /**
      * Destructor.
      */
      ~System();

      /**
      * Read input parameters. 
      */
      void readParameters(std::istream& in);

      /**
      * Get Mixture by reference.
      */
      Mixture& mixture();
      
      /**
      * Get Iterator by reference.
      */
      Iterator& iterator();

   private:

      // Associated Mixture object.
      Mixture mixture_;

      // Pointer to associated iterator.
      Iterator* iteratorPtr_;
   };

   // Inline functions

   inline
   Mixture& System::mixture()
   {
      UTIL_ASSERT(mixturePtr_);
      return mixture_;
   }

   inline
   Iterator& System::iterator()
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }

} // namespace Fd1d
} // namespace Pscf
#endif
