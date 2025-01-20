#ifndef RPG_EXPLICIT_BD_STEP_TPP
#define RPG_EXPLICIT_BD_STEP_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ExplicitBdStep.h"

#include <rpg/fts/brownian/BdSimulator.h>
#include <rpg/fts/compressor/Compressor.h>
#include <rpg/System.h>
#include <pscf/math/IntVec.h>
#include <pscf/cuda/VecOp.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

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
      IntVec<D> dimensions = system().domain().mesh().dimensions();
      int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(w_.capacity() == nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         UTIL_CHECK(w_[i].capacity() == meshSize);
      }
      UTIL_CHECK(dwc_.capacity() == meshSize);
      gaussianField_.allocate(dimensions);
   }

   template <int D>
   bool ExplicitBdStep<D>::step()
   {
      // Array sizes and indices
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int i, j;
      
      // Save current state
      simulator().saveState();

      // Copy current W fields from parent system into w_
      for (i = 0; i < nMonomer; ++i) {
         VecOp::eqV(w_[i], system().w().rgrid(i));
      }
      
      // Constants for dynamics
      const double vSystem = system().domain().unitCell().volume();
      double a = -1.0*mobility_;
      double b = sqrt(2.0*mobility_*double(meshSize)/vSystem);
      
      // Constants for normal distribution
      double stddev = 1.0;
      double mean = 0;
      
      // Modify local field copy wc_
      // Loop over eigenvectors of projected chi matrix
      double evec;
      for (j = 0; j < nMonomer - 1; ++j) {
         RField<D> const & dc = simulator().dc(j);
         
         // Generate normal distributed random floating point numbers
         cudaRandom().normal(gaussianField_, stddev, mean);
         
         // dwc
         VecOp::addVcVc(dwc_, dc, a, gaussianField_, b);
         
         // Loop over monomer types
         for (i = 0; i < nMonomer; ++i) {
            RField<D> & w = w_[i];
            evec = simulator().chiEvecs(j,i);
            VecOp::addEqVc(w, dwc_, evec);
         }
         
      }

      // Set modified fields in parent system
      system().setWRGrid(w_);
      simulator().clearData();

      // Enforce incompressibility (also solves MDE repeatedly)
      bool isConverged = false;
      int compress = simulator().compressor().compress();
      if (compress != 0){
         simulator().restoreState();
      } else {
         isConverged = true;
         UTIL_CHECK(system().hasCFields());
         
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
