/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "MaxOrderParameter.h"

#include <rpg/system/System.h>
#include <rpg/fts/simulator/Simulator.h>
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/field/WFields.h>
#include <rpg/field/CFields.h>

#include <prdc/field/cuda/FFT.h>
#include <prdc/field/cuda/RField.h>
#include <pscf/cuda/VecOp.h>

#include <rp/fts/analyzer/MaxOrderParameterBase.tpp>

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;

   /*
   * Constructor.
   */
   template <int D>
   MaxOrderParameter<D,CUT>::MaxOrderParameter(
                           Simulator<D,CUT>& simulator,
                           System<D,CUT>& system)
    : MaxOrderParameterBase<D,CUT>(simulator, system)
   {}

   /*
   * Compute and return maximum of square magnitude Fourier amplitude.
   */
   template <int D>
   void MaxOrderParameter<D,CUT>::setup()
   {
      // Setup base class
      Base::setup();

      // Allocate psiHost_ array 
      int kSize = Base::kSize_;
      UTIL_CHECK(kSize > 0);
      if (!psiHost_.isAllocated()) {
         psiHost_.allocate(kSize);
      }
      UTIL_CHECK(psiHost_.capacity() == kSize);
   }
      

   /*
   * Compute and return maximum of square magnitude Fourier amplitude.
   */
   template <int D>
   double MaxOrderParameter<D,CUT>::compute()
   {
      // Compute device array psi_ of squared Fourier magnitudes
      Base::computePsi();

      // Copy device array psi_ to corresponing host array psiHost_
      psiHost_ = Base::psi_;

      // Compute maximum from host array
      Base::findMaximum(psiHost_);

      return Base::maxPsi_;
   }

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class MaxOrderParameterBase<1,CUT>;
      template class MaxOrderParameterBase<2,CUT>;
      template class MaxOrderParameterBase<3,CUT>;
      template class MaxOrderParameter<1,CUT>;
      template class MaxOrderParameter<2,CUT>;
      template class MaxOrderParameter<3,CUT>;
   }
}
