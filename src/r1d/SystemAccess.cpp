/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SystemAccess.h"                        // header

namespace Pscf {
namespace R1d {

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

} // namespace R1d
} // namespace Pscf
