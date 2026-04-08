/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MaxOrderParameter.h"

#include <rpc/system/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <prdc/cpu/FFT.h>
#include <prdc/cpu/RField.h>
#include <pscf/cpu/VecOpCx.h>

#include <rp/fts/analyzer/MaxOrderParameter.tpp>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D>
   MaxOrderParameter<D>::MaxOrderParameter(Simulator<D>& simulator,
                                           System<D>& system)
    : Rp::MaxOrderParameter<D, Types<D> >(simulator, system)
   {}

   /*
   * Compute and return maximum of square magnitude Fourier amplitude.
   */
   template <int D>
   double MaxOrderParameter<D>::compute()
   {
      RpMaxOrderParameter::computePsi();
      RpMaxOrderParameter::findMaximum(RpMaxOrderParameter::psi_);
      return RpMaxOrderParameter::maxPsi_;
   }

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class MaxOrderParameter<1, Rpc::Types<1> >;
      template class MaxOrderParameter<2, Rpc::Types<2> >;
      template class MaxOrderParameter<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class MaxOrderParameter<1>;
      template class MaxOrderParameter<2>;
      template class MaxOrderParameter<3>;
   }
}
