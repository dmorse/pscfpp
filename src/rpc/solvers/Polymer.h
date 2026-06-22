#ifndef RPC_POLYMER_H
#define RPC_POLYMER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/Polymer.h>    // base class template
#include <rpc/system/Types.h>      // base class template argument

// Forward declarations
namespace Pscf {
   namespace Rp {
      template <int D, class T> class Propagator;
   }
   namespace Rpc {
      template <int D> class Block;
   }
}

// Explicit instantiation declarations
namespace Pscf {
   extern template class PolymerTmpl< Rpc::Block<1>, Rp::Propagator<1, Rpc::Types<1> > >;
   extern template class PolymerTmpl< Rpc::Block<2>, Rp::Propagator<2, Rpc::Types<2> > >;
   extern template class PolymerTmpl< Rpc::Block<3>, Rp::Propagator<3, Rpc::Types<3> > >;
   namespace Rp {
      extern template class Polymer<1, Rpc::Types<1> >;
      extern template class Polymer<2, Rpc::Types<2> >;
      extern template class Polymer<3, Rpc::Types<3> >;
   }
}
#endif
