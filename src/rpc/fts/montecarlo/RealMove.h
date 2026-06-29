#ifndef RPC_REAL_MOVE_H
#define RPC_REAL_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/RealMove.h>    // base class template
#include <rpc/system/Types.h>              // template argument
#include <rpc/fts/montecarlo/McMove.h>     // indirect base class
#include <prdc/field/cpu/RField.h>               // base class member

namespace Pscf {
namespace Rpc {

   template <int D> class McSimulator;

   using namespace Util;
   using namespace Prdc;

   /**
   * RealMove generates spatially uncorrelated random field changes.
   *
   * Specializations of this template with D = 1, 2, and 3 are derived from 
   * specializations of the base class template Rp::RealMove, and inherit 
   * their public interface and almost all of their source code from this 
   * base class.  See the documentation of this base class template for 
   * details. 
   *
   * \see \ref rp_RealMove_page "Manual Page"
   * \ingroup Rpc_Fts_MonteCarlo_Module
   */
   template <int D>
   class RealMove : public Rp::RealMove<D, Types<D> >
   {
   public:

      /**
      * Constructor.
      *
      * \param simulator  parent McSimulator object
      */
      RealMove(Rp::McSimulator<D, Rpc::Types<D> >& simulator);

   };

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class RealMove<1, Rpc::Types<1> >;
      extern template class RealMove<2, Rpc::Types<2> >;
      extern template class RealMove<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class RealMove<1>;
      extern template class RealMove<2>;
      extern template class RealMove<3>;
   }
}
#endif
