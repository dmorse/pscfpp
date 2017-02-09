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
   * Constructor.
   */
   SystemAccess::SystemAccess(System& system)
    : system_(system)
   {}

   /**
   * Destructor.
   */
   SystemAccess::~SystemAccess()
   {}

} // namespace Fd1d
} // namespace Pscf
