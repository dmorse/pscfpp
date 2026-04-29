/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ConcentrationWriter.h"

#include <rpc/fts/simulator/Simulator.h>
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/FieldIo.h>
#include <rpc/field/CFields.h>
#include <rpc/field/WFields.h>
#include <prdc/cpu/RField.h>

#include <rp/fts/analyzer/ConcentrationWriter.tpp>

namespace Pscf {
namespace Rpc {

   // Constructor.
   template <int D>
   ConcentrationWriter<D>::ConcentrationWriter(
                                  Simulator<D>& simulator,
                                  System<D>& system)
    : Rp::ConcentrationWriter< D, Types<D> >(simulator, system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ConcentrationWriter<1, Rpc::Types<1> >;
      template class ConcentrationWriter<2, Rpc::Types<2> >;
      template class ConcentrationWriter<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class ConcentrationWriter<1>;
      template class ConcentrationWriter<2>;
      template class ConcentrationWriter<3>;
   }
}
