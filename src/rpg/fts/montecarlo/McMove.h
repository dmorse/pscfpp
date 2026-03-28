#ifndef RPG_MC_MOVE_H
#define RPG_MC_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/montecarlo/McMove.h>
#include <rpg/system/Types.h>

namespace Pscf {
namespace Rpg {

   /**
   * McMove is an abstract base class for Monte Carlo moves.
   *
   * The virtual move() function must generate a trial move, decide whether
   * to accept or reject it, and update the associated System fields if
   * it is accepted.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of base class template Rp::McMove, and
   * inherit their public interface and almost all of their source code
   * from this base class.  
   *
   * \see Rp::McMove
   * \ingroup Rpg_Fts_MonteCarlo_Module
   */
   template <int D>
   class McMove : public Rp::McMove<D, Types<D> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Rpg::McSimulator<D> object
      */
      McMove(McSimulator<D>& simulator);

      McMove() = delete;
      McMove(McMove<D> const &) = delete;

   };


}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class McMove<1, Rpg::Types<1> >;
      extern template class McMove<2, Rpg::Types<2> >;
      extern template class McMove<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class McMove<1>;
      extern template class McMove<2>;
      extern template class McMove<3>;
   }
}
#endif
