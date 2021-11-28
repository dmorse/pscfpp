#ifndef PSPC_LINEAR_SWEEP_TPP
#define PSPC_LINEAR_SWEEP_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pspc/System.h>
#include "LinearSweep.h"

namespace Pscf {
namespace Pspc {

   using namespace Util;

   template <int D>
   LinearSweep<D>::LinearSweep(System<D>& system)
   : Sweep<D>(system)
   {}

   template <int D>
   void LinearSweep<D>::readParameters(std::istream& in)
   {
      SweepTmpl< BasisFieldState<D> >::readParameters(in);
      this->read(in, "nParameter", nParameter_);
      parameters_.allocate(nParameter_);
      this->template readDArray< LinearSweepParameter<D> >(in, "parameters", parameters_, nParameter_);
      // verify net zero change in volume fractions
      double sum = 0;
      for (int i = 0; i < nParameter_; ++i) {
         if (parameters_[i].type() == "phi") {
            sum += parameters_[i].change();
         }
      }
      UTIL_CHECK(sum > -0.0001);
      UTIL_CHECK(sum < 0.0001);
   }

   template <int D>
   void LinearSweep<D>::setup()
   {
      for (int i = 0; i < nParameter_; ++i) {
         parameters_[i].setSystem(Sweep<D>::system());
         parameters_[i].getInitial();
      }
   }

   template <int D>
   void LinearSweep<D>::setParameters(double s)
   {
      for (int i = 0; i < nParameter_; ++i) {
         parameters_[i].update(s);
      }
   }

   template <int D>
   void LinearSweep<D>::outputSummary(std::ostream& out)
   {}

}
}

#endif