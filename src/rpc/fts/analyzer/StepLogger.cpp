/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "StepLogger.h"
#include <rp/fts/analyzer/StepLogger.tpp>

namespace Pscf {
namespace Rpc {

   // Constructor.
   template <int D>
   StepLogger<D>::StepLogger(Rp::Simulator<D, Rpc::Types<D> >& simulator, Rp::System<D, Rpc::Types<D> >& system)
    : Rp::StepLogger< D, Types<D> >(simulator, system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class StepLogger< 1, Rpc::Types<1> >;
      template class StepLogger< 2, Rpc::Types<2> >;
      template class StepLogger< 3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class StepLogger<1>;
      template class StepLogger<2>;
      template class StepLogger<3>;
   }
}
