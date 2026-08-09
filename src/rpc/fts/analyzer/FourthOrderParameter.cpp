/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FourthOrderParameter.h"

#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>

#include <prdc/field/cpu/FFT.h>

#include <pscf/cpu/VecOpCx.h>
#include <pscf/cpu/ReduceCx.h>

#include <rp/fts/analyzer/FourthOrderParameterBase.tpp>

namespace Pscf {
namespace Rp {

   /*
   * Constructor.
   */
   template <int D>
   FourthOrderParameter<D,CPT>::FourthOrderParameter(
                                   Simulator<D,CPT>& simulator,
                                   System<D,CPT>& system)
    : Base(simulator, system)
   {}

   /*
   * Initialize Base::prefactor_ protected member variable.
   */
   template <int D>
   void FourthOrderParameter<D,CPT>::computePrefactor()
   {  Base::computePrefactor(Base::prefactor_); }

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class FourthOrderParameterBase<1,CPT>;
      template class FourthOrderParameterBase<2,CPT>;
      template class FourthOrderParameterBase<3,CPT>;
      template class FourthOrderParameter<1,CPT>;
      template class FourthOrderParameter<2,CPT>;
      template class FourthOrderParameter<3,CPT>;
   }
}
