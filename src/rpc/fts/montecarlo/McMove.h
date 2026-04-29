#ifndef RPC_MC_MOVE_H
#define RPC_MC_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/McMove.h>
#include <rpc/system/Types.h>

namespace Pscf {
namespace Rpc {

   /**
   * McMove is an abstract base class for Monte Carlo moves.
   *
   * The virtual move() function must generate a trial move, decide 
   * whether to accept or reject it, and update the associated System 
   * fields if it is accepted.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * specializations of the base class template Rp::McMove< D, Types<D> >, 
   * and inherit their public interface and almost all of their source 
   * code from this base class.  
   *
   * \see Rp::McMove
   * \see \ref psfts_algorithm_MC_sec "Manual Page"
   * \ingroup Rpc_Fts_MonteCarlo_Module
   */
   template <int D>
   class McMove : public Rp::McMove<D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Rpc::McSimulator<D> object
      */
      McMove(McSimulator<D>& simulator);

      /**
      * Destructor
      */
      virtual ~McMove()
      {}

      // Not copyable.
      McMove(McMove<D> const &) = delete;

   };


}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class McMove<1, Rpc::Types<1> >;
      extern template class McMove<2, Rpc::Types<2> >;
      extern template class McMove<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      extern template class McMove<1>;
      extern template class McMove<2>;
      extern template class McMove<3>;
   }
}
#endif
