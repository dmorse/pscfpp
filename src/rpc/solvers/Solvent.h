#ifndef RPC_SOLVENT_H
#define RPC_SOLVENT_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/Solvent.h>    // base class template
#include <rpc/system/Types.h>      // base class template parameter
#include <prdc/cpu/RField.h>       // member of base class

namespace Pscf {
namespace Rpc {

   /**
   * Solver and descriptor for a solvent species.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of base class template Rp::Solvent, 
   * and inherit their public interface and all of their source code from
   * from this base class.  
   *
   * \see Rp::Solvent
   * \ref user_param_solvent_sec "Manual Page"
   * \ingroup Rpc_Solver_Module
   */
   template <int D>
   class Solvent : public Rp::Solvent<D, Rpc::Types<D> >
   {};

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Solvent<1, Rpc::Types<1> >;
      extern template class Solvent<2, Rpc::Types<2> >;
      extern template class Solvent<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class Solvent<1>;
      extern template class Solvent<2>;
      extern template class Solvent<3>;
   }
}
#endif
