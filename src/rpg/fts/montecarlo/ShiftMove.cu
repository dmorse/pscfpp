/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ShiftMove.h"

#include <rpg/fts/montecarlo/McSimulator.h>
#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/WFields.h>
#include <pscf/mesh/Mesh.h>

#include <rp/fts/montecarlo/ShiftMoveBase.tpp>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /*
   * Constructor.
   */
   template <int D>
   ShiftMove<D, Rpg::Types<D> >::ShiftMove(
                                McSimulator<D, Rpg::Types<D> >& simulator)
    : ShiftMoveBaseT(simulator)
   {}

   /*
   * Setup immediately before starting a simulation.
   */
   template <int D>
   void ShiftMove<D, Rpg::Types<D> >::setup()
   {
      ShiftMoveBaseT::setup();

      // Allocate CPU work space
      if (!wOld_.isAllocated()) {
	 UTIL_CHECK(!wNew_.isAllocated()); 
         const int meshSize = system().domain().mesh().size();
         wOld_.allocate(meshSize);
         wNew_.allocate(meshSize);
      }
   }
   /*
   * Compute and store array w_ of shifted fields.
   */
   template <int D>
   void ShiftMove<D, Rpg::Types<D> >::shiftFields(IntVec<D> const & shift)
   {
      IntVec<D> const& dimensions = system().domain().mesh().dimensions();
      const int nMonomer = system().mixture().nMonomer();

      for (int j = 0; j< nMonomer; ++j) {
         wOld_ = system().w().rgrid(j);
	 ShiftMoveBaseT::shiftField(wNew_, wOld_, shift, dimensions);
         ShiftMoveBaseT::w_[j] = wNew_;
      }
   }

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class ShiftMoveBase<1, Rpg::Types<1> >;
      template class ShiftMoveBase<2, Rpg::Types<2> >;
      template class ShiftMoveBase<3, Rpg::Types<3> >;
      template class ShiftMove<1, Rpg::Types<1> >;
      template class ShiftMove<2, Rpg::Types<2> >;
      template class ShiftMove<3, Rpg::Types<3> >;
   }
}
