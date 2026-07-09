/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "IntraCorrelation.h"

#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>

#include <prdc/field/cpu/FFT.h>
#include <prdc/field/cpu/RField.h>
#include <prdc/field/cpu/send.h>

#include <rp/fts/compressor/IntraCorrelation.tpp>

#if 0
namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;

   /*
   * Constructor.
   */
   template <int D>
   Rp::IntraCorrelation<D, Rpc::Types<D> >::IntraCorrelation(
              Rp::System<D, Rpc::Types<D> > const & system)
    : Rp::IntraCorrelation<D, Types<D> >(system)
   {}

}
}
#endif

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class IntraCorrelation<1, Rpc::Types<1> >;
      template class IntraCorrelation<2, Rpc::Types<2> >;
      template class IntraCorrelation<3, Rpc::Types<3> >;
   }
}
