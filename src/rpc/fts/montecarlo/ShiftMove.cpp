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
#include <pscf/mesh/Mesh.h>

#include <rp/fts/montecarlo/ShiftMoveBase.tpp>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc::Cpu;

   /*
   * Constructor.
   */
   template <int D>
   ShiftMove<D, Rpc::Types<D> >::ShiftMove(
		            McSimulator<D, Rpc::Types<D> >& simulator)
    : ShiftMoveBaseT(simulator)
   {}

   /*
   * Compute and store array w_ of shifted fields.
   */
   template <int D>
   void ShiftMove<D, Rpc::Types<D> >::shiftFields(IntVec<D> const & shift)
   {
      IntVec<D> const& dimensions = system().domain().mesh().dimensions();
      const int nMonomer = system().mixture().nMonomer();
      for (int j = 0; j< nMonomer; ++j) {
         RField<D> const & wOld = system().w().rgrid(j);
         RField<D> & wNew = ShiftMoveBaseT::w_[j];
	 ShiftMoveBaseT::shiftField(wNew, wOld, shift, dimensions);
      }
   }

}
}

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class ShiftMoveBase<1, Rpc::Types<1> >;
      template class ShiftMoveBase<2, Rpc::Types<2> >;
      template class ShiftMoveBase<3, Rpc::Types<3> >;
      template class ShiftMove<1, Rpc::Types<1> >;
      template class ShiftMove<2, Rpc::Types<2> >;
      template class ShiftMove<3, Rpc::Types<3> >;
   }
}
