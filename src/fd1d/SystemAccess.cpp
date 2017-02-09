/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SystemAccess.h"                        // header

namespace Pscf {
namespace Fd1d {

   using namespace Util;

   /*
   * Default constructor.
   */
   SystemAccess::SystemAccess()
    : systemPtr_(0)
   {}

   /*
   * Constructor.
   */
   SystemAccess::SystemAccess(System& system)
    : systemPtr_(&system)
   {}

   /**
   * Destructor.
   */
   SystemAccess::~SystemAccess()
   {}

   /*
   * Set the system pointer.
   */
   void SystemAccess::setSystem(System& system)
   {  systemPtr_ = &system; }

} // namespace Fd1d
} // namespace Pscf
