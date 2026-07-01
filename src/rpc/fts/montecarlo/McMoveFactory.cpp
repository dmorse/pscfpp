/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMoveFactory.h"
//#include <rpc/fts/montecarlo/McSimulator.h>

// Subclasses of McMove
#include "RealMove.h"
#include "ForceBiasMove.h"
#include "BdMove.h"
#include "ShiftMove.h"

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   McMoveFactory<D>::McMoveFactory(Rp::McSimulator<D, Rpc::Types<D> >& simulator)
    : simulatorPtr_(&simulator)
   {}

   /*
   * Return a pointer to a instance of McMove subclass className.
   */
   template <int D>
   Rp::McMove<D, Rpc::Types<D> >* McMoveFactory<D>::factory(const std::string &className) const
   {
      Rp::McMove<D, Rpc::Types<D> >* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      // Try to match classname
      if (className == "RealMove") {
         ptr = new Rp::RealMove<D, Rpc::Types<D> >(*simulatorPtr_);
      } else if (className == "ForceBiasMove") {
         ptr = new Rp::ForceBiasMove<D, Rpc::Types<D> >(*simulatorPtr_);
      } else if (className == "BdMove") {
         ptr = new Rp::BdMove<D, Rpc::Types<D> >(*simulatorPtr_);
      } else if (className == "ShiftMove") {
         ptr = new Rp::ShiftMove<D, Rpc::Types<D> >(*simulatorPtr_);
      }

      return ptr;
   }

   // Explicit instantiation definitions
   template class McMoveFactory<1>;
   template class McMoveFactory<2>;
   template class McMoveFactory<3>;
}
}
