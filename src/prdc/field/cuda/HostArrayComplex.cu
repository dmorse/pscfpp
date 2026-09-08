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
   * Default constructor.
   */
   HostArrayComplex::HostArrayComplex() 
    : HostArray<cudaComplex,CUT>()
   {}

   /*
   * Allocating constructor.
   */
   HostArrayComplex::HostArrayComplex(int capacity)
    : HostArray<cudaComplex,CUT>(capacity)
   {}

   /*
   * Copy constructor.
   */
   HostArrayComplex::HostArrayComplex(HostArrayComplex const& other)
    : HostArray<cudaComplex,CUT>(other)
   {}

   /*
   * Destructor.
   */
   HostArrayComplex::~HostArrayComplex()
   {}

} // namespace Prdc
} // namespace Pscf
