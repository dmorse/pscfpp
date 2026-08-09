/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ShiftMove_u.h"
#include <rp/fts/montecarlo/ShiftMoveBase.tpp>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /*
   * Constructor.
   */
   template <int D>
   ShiftMove<D,CUT>::ShiftMove(
                                McSimulator<D,CUT>& simulator)
    : ShiftMoveBaseT(simulator)
   {}

   /*
   * Setup immediately before starting a simulation.
   */
   template <int D>
   void ShiftMove<D,CUT>::setup()
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
   void ShiftMove<D,CUT>::shiftFields(IntVec<D> const & shift)
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
      template class ShiftMoveBase<1,CUT>;
      template class ShiftMoveBase<2,CUT>;
      template class ShiftMoveBase<3,CUT>;
      template class ShiftMove<1,CUT>;
      template class ShiftMove<2,CUT>;
      template class ShiftMove<3,CUT>;
   }
}
