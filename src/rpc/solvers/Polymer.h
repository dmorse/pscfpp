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

namespace Pscf {
namespace Rpc {

   /**
   * Descriptor and solver for one polymer species.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of base class template Rp::Polymer, and
   * inherit their public interface and all of their source code from this
   * base class.  
   *
   * \see Rp::Polymer
   * \ref user_param_polymer_sec "Manual Page"
   * \ingroup Rpc_Solver_Module
   */
   template <int D>
   class Polymer : public Rp::Polymer< D, Types<D> >
   {};

}
}

// Explicit instantiation declarations
namespace Pscf {
   extern template class PolymerTmpl< Rpc::Block<1>, Rpc::Propagator<1> >;
   extern template class PolymerTmpl< Rpc::Block<2>, Rpc::Propagator<2> >;
   extern template class PolymerTmpl< Rpc::Block<3>, Rpc::Propagator<3> >;
   namespace Rp {
      extern template class Polymer<1, Rpc::Types<1> >;
      extern template class Polymer<2, Rpc::Types<2> >;
      extern template class Polymer<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class Polymer<1>;
      extern template class Polymer<2>;
      extern template class Polymer<3>;
   }
}
#endif
