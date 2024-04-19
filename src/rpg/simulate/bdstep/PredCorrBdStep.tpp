#ifndef RPG_PRED_CORR_BD_STEP_TPP
#define RPG_PRED_CORR_BD_STEP_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PredCorrBdStep.h"

#include <rpg/simulate/bdstep/BdSimulator.h>
#include <rpg/simulate/compressor/Compressor.h>
#include <rpg/System.h>
#include <pscf/math/IntVec.h>
#include <util/random/Random.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /*
   * Constructor.
   */
   template <int D>
   PredCorrBdStep<D>::PredCorrBdStep(BdSimulator<D>& simulator)
    : BdStep<D>(simulator),
      wp_(),
      wf_(),
      dci_(),
      eta_(),
      dwc_(),
      dwp_(),
      mobility_(0.0)
   {}

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   PredCorrBdStep<D>::~PredCorrBdStep()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void PredCorrBdStep<D>::readParameters(std::istream &in)
   {
      read(in, "mobility", mobility_);

      // Allocate memory for private containers
      int nMonomer = system().mixture().nMonomer();
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      wp_.allocate(nMonomer);
      wf_.allocate(nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         wp_[i].allocate(meshDimensions);
         wf_[i].allocate(meshDimensions);
      }
      dci_.allocate(nMonomer-1);
      eta_.allocate(nMonomer-1);
      for (int i=0; i < nMonomer - 1; ++i) {
         dci_[i].allocate(meshDimensions);
         eta_[i].allocate(meshDimensions);
      }
      dwc_.allocate(meshDimensions);
      dwp_.allocate(meshDimensions);
   }

   template <int D>
   void PredCorrBdStep<D>::setup()
   {
      // Check array capacities
      int meshSize = system().domain().mesh().size();
      int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(wp_.capacity() == nMonomer);
      UTIL_CHECK(wf_.capacity() == nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         UTIL_CHECK(wp_[i].capacity() == meshSize);
         UTIL_CHECK(wf_[i].capacity() == meshSize);
      }
      UTIL_CHECK(dci_.capacity() == nMonomer-1);
      UTIL_CHECK(eta_.capacity() == nMonomer-1);
      for (int i=0; i < nMonomer - 1; ++i) {
         UTIL_CHECK(dci_[i].capacity() == meshSize);
         UTIL_CHECK(eta_[i].capacity() == meshSize);
      }
      UTIL_CHECK(dwc_.capacity() == meshSize);
      UTIL_CHECK(dwp_.capacity() == meshSize);
   }

   template <int D>
   bool PredCorrBdStep<D>::step()
   {
      // Array sizes and indices
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int i, j;
      
      // Save current state
      simulator().saveState();
      
      // GPU resources  
      int nBlocks, nThreads;
      ThreadGrid::setThreadsLogical(meshSize, nBlocks, nThreads);
      
      // Copy current W fields from parent system
      DArray<RField<D>> const * Wr = &system().w().rgrid();
      for (i = 0; i < nMonomer; ++i) {
         assignReal<<<nBlocks, nThreads>>>
            (wp_[i].cField(), (*Wr)[i].cField(), meshSize);
         assignReal<<<nBlocks, nThreads>>>
            (wf_[i].cField(), wp_[i].cField(), meshSize);
      }

      // Store initial value of pressure field
      assignReal<<<nBlocks, nThreads>>>
         (dwp_.cField(), simulator().wc(nMonomer-1).cField(), meshSize);

      // Constants for dynamics
      const double vSystem = system().domain().unitCell().volume();
      double a = -1.0*mobility_;
      double b = sqrt(2.0*mobility_*double(meshSize)/vSystem);
      
      // Constants for normal distribution
      double stddev = 1.0;
      double mean = 0;
      
      // Construct all random displacement components
      for (j = 0; j < nMonomer - 1; ++j) {
         RField<D> & eta = eta_[j];
         cudaRandom().normal(eta.cField(), meshSize, (cudaReal)stddev, (cudaReal)mean);
         scaleReal<<<nBlocks, nThreads>>>(eta.cField(), b, meshSize);
      }

      // Compute predicted state wp_, and store initial force dci_

      // Loop over eigenvectors of projected chi matrix
      double evec;
      for (j = 0; j < nMonomer - 1; ++j) {
         RField<D> const & dc = simulator().dc(j);
         RField<D> const & eta = eta_[j];
         RField<D> & dci = dci_[j];
         
         // dwc_[k] = a*dc[k] + eta[k];
         assignReal<<<nBlocks, nThreads>>>(dwc_.cField(), eta.cField(), meshSize);
         pointWiseAddScale<<<nBlocks, nThreads>>>(dwc_.cField(), dc.cField(), a, meshSize);
         
         // dci[k] = dc[k];
         assignReal<<<nBlocks, nThreads>>>(dci.cField(), dc.cField(), meshSize);

         // Loop over monomer types
         for (i = 0; i < nMonomer; ++i) {
            RField<D> & wp = wp_[i];
            evec = simulator().chiEvecs(j,i);
            pointWiseAddScale<<<nBlocks, nThreads>>>(wp.cField(), dwc_.cField(), evec, meshSize);
         }
      }

      // Set modified fields at predicted state wp_
      system().setWRGrid(wp_);
      
      // Enforce incompressibility (also solves MDE repeatedly)
      bool isConverged = false;
      int compress = simulator().compressor().compress();
      if (compress != 0){
         simulator().restoreState();
      } else {
         UTIL_CHECK(system().hasCFields());

         // Compute components and derivatives at wp_
         simulator().clearData();
         simulator().computeWc();
         simulator().computeCc();
         simulator().computeDc();
         
         // Compute change in pressure field
         RField<D> const & wp = simulator().wc(nMonomer-1);
         
         // dwp_[k] = wp[k] - dwp_[k]
         pointWiseBinarySubtract<<<nBlocks, nThreads>>>(wp.cField(), dwp_.cField(), dwp_.cField(), meshSize);

         // Adjust pressure field
         for (i = 0; i < nMonomer; ++i) {
            RField<D> & wf = wf_[i];
            pointWiseAdd<<<nBlocks, nThreads>>>(wf.cField(), dwp_.cField(), meshSize);
         }
         
         // Full step (corrector)
         double ha = 0.5*a;
         for (j = 0; j < nMonomer - 1; ++j) {
            RField<D> const & dcp = simulator().dc(j);
            RField<D> const & dci = dci_[j];
            RField<D> const & eta = eta_[j];
            
            // dwc_[k] = ha*( dci[k] + dcp[k]) + eta[k];
            pointWiseBinaryAdd<<<nBlocks, nThreads>>>(dci.cField(), dcp.cField(), dwc_.cField(), meshSize);
            scaleReal<<<nBlocks, nThreads>>>(dwc_.cField(), ha, meshSize);
            pointWiseAdd<<<nBlocks, nThreads>>>(dwc_.cField(), eta.cField(), meshSize);

            for (i = 0; i < nMonomer; ++i) {
               RField<D> & wf = wf_[i];
               evec = simulator().chiEvecs(j,i);
               pointWiseAddScale<<<nBlocks, nThreads>>>(wf.cField(), dwc_.cField(), evec, meshSize);
            }
         }
         
         // Set fields at final point
         system().setWRGrid(wf_);
         
         // Enforce incompressibility for final point
         int compress2 = simulator().compressor().compress();
         if (compress2 != 0){
            simulator().restoreState();
         }  else {
            isConverged = true;
            UTIL_CHECK(system().hasCFields());

            // Compute components and derivatives at final point
            simulator().clearState();
            simulator().clearData();
            simulator().computeWc();
            simulator().computeCc();
            simulator().computeDc();
            
         }
      }
      
      return isConverged;
   }

}
}
#endif
