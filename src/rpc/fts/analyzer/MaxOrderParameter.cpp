/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MaxOrderParameter.h"

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>

#include <prdc/field/cpu/FFT.h>
#include <prdc/field/cpu/RField.h>
#include <pscf/cpu/VecOpCx.h>

#include <rp/fts/analyzer/MaxOrderParameterBase.tpp>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D>
   MaxOrderParameter<D,CPT>::MaxOrderParameter(
                           Simulator<D,CPT>& simulator,
                           System<D,CPT>& system)
    : MaxOrderParameterBase<D,CPT>(simulator, system)
   {}

   /*
   * Compute and return maximum of square magnitude Fourier amplitude.
   */
   template <int D>
   double MaxOrderParameter<D,CPT>::compute()
   {
      Base::computePsi();
      Base::findMaximum(Base::psi_);
      return Base::maxPsi_;
   }

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class MaxOrderParameterBase<1,CPT>;
      template class MaxOrderParameterBase<2,CPT>;
      template class MaxOrderParameterBase<3,CPT>;
      template class MaxOrderParameter<1,CPT>;
      template class MaxOrderParameter<2,CPT>;
      template class MaxOrderParameter<3,CPT>;
   }
}
