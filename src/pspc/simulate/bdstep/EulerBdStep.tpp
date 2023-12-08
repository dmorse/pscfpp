#ifndef PSPC_EULER_BD_STEP_TPP
#define PSPC_EULER_BD_STEP_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "EulerBdStep.h"

#include <pspc/simulate/BdSimulator.h>
#include <pspc/System.h>
#include <pscf/math/IntVec.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;
   using namespace Prdc::Cpu;

   /*
   * Constructor.
   */
   template <int D>
   EulerBdStep<D>::EulerBdStep(BdSimulator<D>& bdSimulator)
    : BdStep<D>(bdSimulator),
      w_(),
      mobility_(0.0)
   {}

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   EulerBdStep<D>::~EulerBdStep()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void EulerBdStep<D>::readParameters(std::istream &in)
   {
      read<double>(in, "mobility", mobility_);

      IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      int nMonomer = system().mixture().nMonomer();

      w_.allocate(nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         w_[i].allocate(meshDimensions);
      }

   }

   template <int D>
   void EulerBdStep<D>::setup()
   {
       // Check capacity and dimensions of w_;
   }

   template <int D>
   void EulerBdStep<D>::step()
   {
      // Copy W fields from system into w_

      // Add deterministic displacement

      // Generate random numbers
      // Add random displacements

      // Apply new w fields to parent system
   }


}
}
#endif
