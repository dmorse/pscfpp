/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HostArrayComplex.h"

namespace Pscf {
namespace Prdc {

   using namespace Util;

   /*
   * Allocating constructor.
   */
   HostArrayComplex::HostArrayComplex(int capacity)
    : HostArray<cudaComplex,CUT>(capacity)
   {}

} // namespace Prdc
} // namespace Pscf
