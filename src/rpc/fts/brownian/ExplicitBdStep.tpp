#ifndef RPC_EXPLICIT_BD_STEP_TPP
#define RPC_EXPLICIT_BD_STEP_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ExplicitBdStep.h"
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/fts/brownian/BdSimulator.h>
#include <rpc/fts/compressor/Compressor.h>
#include <pscf/math/IntVec.h>
#include <util/random/Random.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Prdc::Cpu;

   /*
   * Constructor.
   */
   template <int D>
   ExplicitBdStep<D>::ExplicitBdStep(BdSimulator<D>& simulator)
    : BdStep<D>(simulator),
      w_(),
      dwc_(),
      mobility_(0.0)
   {}

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   ExplicitBdStep<D>::~ExplicitBdStep()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void ExplicitBdStep<D>::readParameters(std::istream &in)
   {
      read(in, "mobility", mobility_);

      // Allocate memory for private containers
      int nMonomer = system().mixture().nMonomer();
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      w_.allocate(nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         w_[i].allocate(meshDimensions);
      }
      dwc_.allocate(meshDimensions);

   }

   template <int D>
   void ExplicitBdStep<D>::setup()
   {
      // Check array capacities
      int meshSize = system().domain().mesh().size();
      int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(w_.capacity() == nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         UTIL_CHECK(w_[i].capacity() == meshSize);
      }
      UTIL_CHECK(dwc_.capacity() == meshSize);
   }

   template <int D>
   bool ExplicitBdStep<D>::step()
   {
      // Array sizes and indices
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int i, j, k;

      // Save current state
      simulator().saveState();

      // Copy current W fields from parent system into w_
      for (i = 0; i < nMonomer; ++i) {
         w_[i] = system().w().rgrid(i);
      }

      // Constants for dynamics
      const double vSystem = system().domain().unitCell().volume();
      const double a = -1.0*mobility_;
      const double b = sqrt(2.0*mobility_*double(meshSize)/vSystem);

      // Modify local field copy wc_
      // Loop over eigenvectors of projected chi matrix
      double dwd, dwr, evec;
      for (j = 0; j < nMonomer - 1; ++j) {
         RField<D> const & dc = simulator().dc(j);
         for (k = 0; k < meshSize; ++k) {
            dwd = a*dc[k];
            dwr = b*random().gaussian();
            dwc_[k] = dwd + dwr;
         }
         // Loop over monomer types
         for (i = 0; i < nMonomer; ++i) {
            RField<D> & w = w_[i];
            evec = simulator().chiEvecs(j,i);
            for (k = 0; k < meshSize; ++k) {
               w[k] += evec*dwc_[k];
            }
         }
      }

      // Set modified fields in parent system
      system().w().setRGrid(w_);
      simulator().clearData();

      // Enforce incompressibility (also solves MDE repeatedly)
      bool isConverged = false;
      int compress = simulator().compressor().compress();
      if (compress != 0){
         simulator().restoreState();
      } else {
         isConverged = true;
         UTIL_CHECK(system().c().hasData());

         // Evaluate component properties in new state
         simulator().clearState();
         simulator().computeWc();
         simulator().computeCc();
         simulator().computeDc();
      }

      return isConverged;
   }

}
}
#endif
