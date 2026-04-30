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

#include <prdc/cuda/FFT.h>
#include <prdc/cuda/RField.h>
#include <pscf/cuda/VecOp.h>

#include <rp/fts/analyzer/MaxOrderParameter.tpp>

namespace Pscf {
namespace Rpg {

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
   void MaxOrderParameter<D>::setup()
   {
      // Setup base class
      RpMaxOrderParameter::setup();

      // Allocate psiHost_ array 
      int kSize = RpMaxOrderParameter::kSize_;
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
   double MaxOrderParameter<D>::compute()
   {
      // Compute device array psi_ of squared Fourier magnitudes
      RpMaxOrderParameter::computePsi();

      // Copy device array psi_ to corresponing host array psiHost_
      psiHost_ = RpMaxOrderParameter::psi_;

      // Compute maximum from host array
      RpMaxOrderParameter::findMaximum(psiHost_);

      return RpMaxOrderParameter::maxPsi_;
   }

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class MaxOrderParameter<1, Rpg::Types<1> >;
      template class MaxOrderParameter<2, Rpg::Types<2> >;
      template class MaxOrderParameter<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      template class MaxOrderParameter<1>;
      template class MaxOrderParameter<2>;
      template class MaxOrderParameter<3>;
   }
}
