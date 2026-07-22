#ifndef RPC_RAMP_H
#define RPC_RAMP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/ramp/Ramp.h>            // base class template
#include <pscf/backends/CPT.h>            // base class template argument

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Ramp<1,CPT>;
      extern template class Ramp<2,CPT>;
      extern template class Ramp<3,CPT>;
   }
}
#endif
