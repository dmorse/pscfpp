/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ShiftMove.h"
#include <rpc/fts/montecarlo/McSimulator.h>
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/WFields.h>

#include <rp/fts/montecarlo/ShiftMove.tpp>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc::Cpu;

   /*
   * Constructor.
   */
   template <int D>
   ShiftMove<D>::ShiftMove(McSimulator<D>& simulator)
    : RpShiftMove(simulator)
   {}

   /*
   * Compute and store array w_ of shifted fields.
   */
   template <int D>
   void ShiftMove<D>::shiftFields(IntVec<D> const & shift)
   {
      IntVec<D> const& dimensions = system().domain().mesh().dimensions();
      const int nMonomer = system().mixture().nMonomer();
      for (int j = 0; j< nMonomer; ++j) {
         RField<D> const & wOld = system().w().rgrid(j);
         RField<D> & wNew = RpShiftMove::w_[j];
	 RpShiftMove::shiftField(wNew, wOld, shift, dimensions);
      }
   }

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class Rp::ShiftMove<1, Rpc::Types<1> >;
      template class Rp::ShiftMove<2, Rpc::Types<2> >;
      template class Rp::ShiftMove<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class ShiftMove<1>;
      template class ShiftMove<2>;
      template class ShiftMove<3>;
   }
}
