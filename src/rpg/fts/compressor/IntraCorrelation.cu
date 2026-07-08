/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IntraCorrelation.h"

#include <rpg/system/System.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>

#include <prdc/field/cuda/FFT.h>
#include <prdc/field/cuda/RField.h>
#include <prdc/field/cuda/send.h>

#include <rp/fts/compressor/IntraCorrelation.tpp>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Prdc;

   /*
   * Constructor.
   */
   template <int D>
   IntraCorrelation<D>::IntraCorrelation(
            Rp::System<D, Rpg::Types<D> > const & system)
    : Rp::IntraCorrelation<D, Types<D> >(system)
   {}

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class IntraCorrelation<1, Rpg::Types<1> >;
      template class IntraCorrelation<2, Rpg::Types<2> >;
      template class IntraCorrelation<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class IntraCorrelation<1>;
      template class IntraCorrelation<2>;
      template class IntraCorrelation<3>;
   }
}
