/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MaxOrderParameter_c.h"
#include <pscf/backend/cpp/VecOpCx.h>

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

   // Explicit instantiation definitions
   template class MaxOrderParameterBase<1,CPT>;
   template class MaxOrderParameterBase<2,CPT>;
   template class MaxOrderParameterBase<3,CPT>;
   template class MaxOrderParameter<1,CPT>;
   template class MaxOrderParameter<2,CPT>;
   template class MaxOrderParameter<3,CPT>;

} // namespace Rp
} // namespace Pscf
