#ifndef PSPC_HAMILTONIAN_AUTO_CORR_TPP
#define PSPC_HAMILTONIAN_AUTO_CORR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HamiltonianAutoCorr.h"

#include <pspc/System.h>
#include <pspc/simulate/Simulator.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   HamiltonianAutoCorr<D>::HamiltonianAutoCorr(Simulator<D>& simulator, 
                                               System<D>& system)
    : Analyzer<D>(),
      simulatorPtr_(&simulator),
      systemPtr_(&(simulator.system())),
      bufferCapacity_(64),
      maxStageId_(0),
      blockFactor_(2),
      isInitialized_(false)
   {  setClassName("HamiltonianAutoCorr"); }

   /*
   * Read interval and outputFileName. 
   */
   template <int D>
   void HamiltonianAutoCorr<D>::readParameters(std::istream& in) 
   {
      // Read interval and parameters for AutoCorrelatio
      readInterval(in);
      readOutputFileName(in);
      // Read parameters for autocorrelation
      readOptional(in,"bufferCapacity", bufferCapacity_);
      readOptional(in,"maxStageId", maxStageId_);
      readOptional(in,"blockFactor", blockFactor_);
   }

   /*
   * HamiltonianAutoCorr setup
   */
   template <int D>
   void HamiltonianAutoCorr<D>::setup() 
   {
      // Set all parameters and allocate to initialize state for autocorrelation 
      accumulator_.setParam(bufferCapacity_, maxStageId_, blockFactor_);
      accumulator_.clear();
      isInitialized_ = true;
   }
   
   /*
   * Sample hamiltonian and add to accumulator.
   */
   template <int D>
   void HamiltonianAutoCorr<D>::sample(long iStep) 
   {
      UTIL_CHECK(simulator().hasWc());
      if (isAtInterval(iStep))  {
         if (!simulator().hasHamiltonian()) {
            simulator().computeHamiltonian();
         }
      double hamiltonian = simulator().hamiltonian();
      accumulator_.sample(hamiltonian);
      }
   }
   
   template <int D>
   void HamiltonianAutoCorr<D>::output()
   {
      system().fileMaster().openOutputFile(outputFileName(), outputFile_);
      accumulator_.output(outputFile_); 
      outputFile_.close();
   } 
   
}
}
#endif
