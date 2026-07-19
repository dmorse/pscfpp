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
#include <rpc/field/WFields.h>

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
   MaxOrderParameter<D, Cpp<D> >::MaxOrderParameter(
                           Simulator<D, Cpp<D> >& simulator,
                           System<D, Cpp<D> >& system)
    : MaxOrderParameterBase<D, Cpp<D> >(simulator, system)
   {}

   /*
   * Compute and return maximum of square magnitude Fourier amplitude.
   */
   template <int D>
   double MaxOrderParameter<D, Cpp<D> >::compute()
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
      template class MaxOrderParameterBase<1, Cpp<1> >;
      template class MaxOrderParameterBase<2, Cpp<2> >;
      template class MaxOrderParameterBase<3, Cpp<3> >;
      template class MaxOrderParameter<1, Cpp<1> >;
      template class MaxOrderParameter<2, Cpp<2> >;
      template class MaxOrderParameter<3, Cpp<3> >;
   }
}
