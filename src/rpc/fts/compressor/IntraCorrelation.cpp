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

#include <prdc/cpu/FFT.h>
#include <prdc/cpu/RField.h>

#include <rp/fts/compressor/IntraCorrelation.tpp>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc;

   /*
   * Constructor.
   */
   template <int D>
   IntraCorrelation<D>::IntraCorrelation(System<D> const & system)
    : Rp::IntraCorrelation<D, Types<D> >(system)
   {}

   #if 0
   /*
   * Compute k-space array of intramolecular correlation functions.
   */
   template<int D>
   void
   IntraCorrelation<D>::computeOmegaTotal(RField<D>& correlations)
   {
      computeMeshProperties();

      // Compute total intramolecular correlation function
      UTIL_CHECK(correlations.capacity() == kSize_);
      if (!correlationMixturePtr_->isAllocated()) {
         correlationMixturePtr_->allocate();
      }
      correlationMixturePtr_->setup();
      correlationMixturePtr_->computeOmegaTotal(Gsq_, correlations);

   }
   #endif

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class IntraCorrelation<1, Rpc::Types<1> >;
      template class IntraCorrelation<2, Rpc::Types<2> >;
      template class IntraCorrelation<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class IntraCorrelation<1>;
      template class IntraCorrelation<2>;
      template class IntraCorrelation<3>;
   }
}
