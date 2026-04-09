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

#include <rp/fts/montecarlo/ShiftMove.tpp>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /*
   * Constructor.
   */
   template <int D>
   ShiftMove<D>::ShiftMove(McSimulator<D>& simulator)
    : RpShiftMove(simulator)
   {}

   /*
   * Setup immediately before starting a simulation.
   */
   template <int D>
   void ShiftMove<D>::setup()
   {
      RpShiftMove::setup();

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
   void ShiftMove<D>::shiftFields(IntVec<D> const & shift)
   {
      IntVec<D> const& dimensions = system().domain().mesh().dimensions();
      const int nMonomer = system().mixture().nMonomer();

      for (int j = 0; j< nMonomer; ++j) {
         wOld_ = system().w().rgrid(j);
	 RpShiftMove::shiftField(wNew_, wOld_, shift, dimensions);
         RpShiftMove::w_[j] = wNew_;
      }
   }

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class Rp::ShiftMove<1, Rpg::Types<1> >;
      template class Rp::ShiftMove<2, Rpg::Types<2> >;
      template class Rp::ShiftMove<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class ShiftMove<1>;
      template class ShiftMove<2>;
      template class ShiftMove<3>;
   }
}
