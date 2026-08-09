/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ShiftMove_c.h"
#include <rp/fts/montecarlo/ShiftMoveBase.tpp>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   ShiftMove<D,CPT>::ShiftMove(
		            McSimulator<D,CPT>& simulator)
    : ShiftMoveBaseT(simulator)
   {}

   /*
   * Compute and store array w_ of shifted fields.
   */
   template <int D>
   void ShiftMove<D,CPT>::shiftFields(IntVec<D> const & shift)
   {
      IntVec<D> const& dimensions = system().domain().mesh().dimensions();
      const int nMonomer = system().mixture().nMonomer();
      for (int j = 0; j< nMonomer; ++j) {
         RField<D,CPT> const & wOld = system().w().rgrid(j);
         RField<D,CPT> & wNew = ShiftMoveBaseT::w_[j];
	 ShiftMoveBaseT::shiftField(wNew, wOld, shift, dimensions);
      }
   }

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class ShiftMoveBase<1,CPT>;
      template class ShiftMoveBase<2,CPT>;
      template class ShiftMoveBase<3,CPT>;
      template class ShiftMove<1,CPT>;
      template class ShiftMove<2,CPT>;
      template class ShiftMove<3,CPT>;
   }
}
