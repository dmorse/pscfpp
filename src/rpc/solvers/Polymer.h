#ifndef RPC_POLYMER_H
#define RPC_POLYMER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/Polymer.h>    // base class template
#include <pscf/backends/CPT.h>      // base class template argument

// Explicit instantiation declarations
namespace Pscf {
   extern template 
   class PolymerTmpl< Rp::Block<1,CPT>, Rp::Propagator<1,CPT> >;
   extern template 
   class PolymerTmpl< Rp::Block<2,CPT>, Rp::Propagator<2,CPT> >;
   extern template 
   class PolymerTmpl< Rp::Block<3,CPT>, Rp::Propagator<3,CPT> >;
   namespace Rp {
      extern template class Polymer<1,CPT>;
      extern template class Polymer<2,CPT>;
      extern template class Polymer<3,CPT>;
   }
}
#endif
