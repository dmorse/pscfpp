/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HostDArrayComplex.h"

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Default constructor.
   */
   HostDArrayComplex::HostDArrayComplex() 
    : HostDArray<cudaComplex>()
   {}

   /*
   * Allocating constructor.
   */
   HostDArrayComplex::HostDArrayComplex(int capacity)
    : HostDArray<cudaComplex>(capacity)
   {}

   /*
   * Copy constructor.
   */
   HostDArrayComplex::HostDArrayComplex(HostDArrayComplex const& other)
    : HostDArray<cudaComplex>(other)
   {}

   /*
   * Destructor.
   */
   HostDArrayComplex::~HostDArrayComplex()
   {}

} // namespace Rpg
} // namespace Pscf
