/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FourthOrderParameter.h"

#include <rpc/system/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>

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
   FourthOrderParameter<D, CppTp<D> >::FourthOrderParameter(
                                   Simulator<D, CppTp<D> >& simulator,
                                   System<D, CppTp<D> >& system)
    : Base(simulator, system)
   {}

   /*
   * Initialize Base::prefactor_ protected member variable.
   */
   template <int D>
   void FourthOrderParameter<D, CppTp<D> >::computePrefactor()
   {  Base::computePrefactor(Base::prefactor_); }

}
}

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class FourthOrderParameterBase<1, CppTp<1> >;
      template class FourthOrderParameterBase<2, CppTp<2> >;
      template class FourthOrderParameterBase<3, CppTp<3> >;
      template class FourthOrderParameter<1, CppTp<1> >;
      template class FourthOrderParameter<2, CppTp<2> >;
      template class FourthOrderParameter<3, CppTp<3> >;
   }
}
