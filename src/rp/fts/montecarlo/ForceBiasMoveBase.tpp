#ifndef RP_FORCE_BIAS_MOVE_BASE_TPP
#define RP_FORCE_BIAS_MOVE_BASE_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ForceBiasMoveBase.h"

#include <rp/fts/montecarlo/McSimulator.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/fts/compressor/Compressor.h>
#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/CFields.h>
#include <rp/field/WFields.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/math/IntVec.h>

#include <util/param/ParamComposite.h>
#include <util/random/Random.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   ForceBiasMoveBase<D,T>::ForceBiasMoveBase(McSimulator<D,T>& simulator)
    : McMove<D,T>(simulator),
      w_(),
      dwc_(),
      mobility_(0.0)
   {  ParamComposite::setClassName("ForceBiasMove"); }

   /*
   * Read body of parameter file block and allocate memory.
   */
   template <int D, class T>
   void ForceBiasMoveBase<D,T>::readParameters(std::istream &in)
   {
      McMove<D,T>::readProbability(in);
      ParamComposite::read(in, "mobility", mobility_);

      // Allocate memory for private containers
      int nMonomer = system().mixture().nMonomer();
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();
      w_.allocate(nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         w_[i].allocate(meshDimensions);
      }
      dc_.allocate(nMonomer-1);
      dwc_.allocate(nMonomer-1);
      for (int i=0; i < nMonomer - 1; ++i) {
         dc_[i].allocate(meshDimensions);
         dwc_[i].allocate(meshDimensions);
      }
      biasField_.allocate(meshDimensions);
      eta_.allocate(meshDimensions);
   }

   /*
   * Setup before entering main simulation loop.
   */
   template <int D, class T>
   void ForceBiasMoveBase<D,T>::setup()
   {
      // Check array capacities
      int meshSize = system().domain().mesh().size();
      int nMonomer = system().mixture().nMonomer();
      UTIL_CHECK(w_.capacity() == nMonomer);
      for (int i=0; i < nMonomer; ++i) {
         UTIL_CHECK(w_[i].capacity() == meshSize);
      }
      UTIL_CHECK(dc_.capacity() == nMonomer-1);
      UTIL_CHECK(dwc_.capacity() == nMonomer-1);
      for (int i = 0; i < nMonomer - 1; ++i) {
         UTIL_CHECK(dc_[i].capacity() == meshSize);
         UTIL_CHECK(dwc_[i].capacity() == meshSize);
      }

      McMove<D,T>::setup();
   }

   /*
   * Attempt and accept or reject MC move
   */
   template <int D, class T>
   bool ForceBiasMoveBase<D,T>::move()
   {
      McMove<D,T>::totalTimer_.start();
      McMove<D,T>::incrementNAttempt();

      // Preconditions
      UTIL_CHECK(simulator().hasWc());
      UTIL_CHECK(simulator().hasDc());
      UTIL_CHECK(simulator().hasHamiltonian());

      // Array sizes and index variables
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int i, j;

      // Get current Hamiltonian
      double oldHamiltonian = simulator().hamiltonian();

      // Save current state
      simulator().saveState();

      // Clear eigen-components of the fields and Hamiltonian
      simulator().clearData();

      McMove<D,T>::attemptMoveTimer_.start();

      // Copy current W fields from parent system into wc_
      for (i = 0; i < nMonomer; ++i) {
         VecOp::eqV(w_[i], system().w().rgrid(i));
      }

      // Copy derivative fields into dc_
      for (i = 0; i < nMonomer - 1; ++i) {
         VecOp::eqV(dc_[i], simulator().dc(i));
      }

      // Constants for dynamics
      const double vSystem = system().domain().unitCell().volume();
      const double vNode = vSystem/double(meshSize);
      const double a = -1.0*mobility_;
      const double b = sqrt(2.0*mobility_/vNode);
      const double stddev = 1.0;
      const double mean = 0.0;

      // Modify local variables dwc_ and w_
      // Loop over composition eigenvectors of projected chi matrix
      for (j = 0; j < nMonomer - 1; ++j) {

         // Generate vector of normal distributed random numbers
         McMove<D,T>::vecRandom().normal(eta_, stddev, mean);

         // Compute vector dwc_[j] of field component changes
         VecOp::addVcVc(dwc_[j], dc_[j], a, eta_, b);

         // Loop over monomer types to add to w_
         for (i = 0; i < nMonomer; ++i) {
            double evec = simulator().chiEvecs(j,i);
            VecOp::addEqVc(w_[i], dwc_[j], evec);
         }

      }

      // Set modified fields in parent system
      system().w().setRGrid(w_);
      simulator().clearData();

      McMove<D,T>::attemptMoveTimer_.stop();

      // Call compressor
      McMove<D,T>::compressorTimer_.start();
      int compress = simulator().compressor().compress();
      UTIL_CHECK(system().c().hasData());
      McMove<D,T>::compressorTimer_.stop();

      bool isConverged = false;
      if (compress != 0){
         McMove<D,T>::incrementNFail();
         simulator().restoreState();
      } else {
         isConverged = true;

         // Compute eigenvector components of current fields
         McMove<D,T>::componentTimer_.start();
         simulator().computeWc();
         UTIL_CHECK(system().c().hasData());
         simulator().computeCc();
         simulator().computeDc();
         McMove<D,T>::componentTimer_.stop();

         // Evaluate new Hamiltonian
         McMove<D,T>::hamiltonianTimer_.start();
         simulator().computeHamiltonian();
         double newHamiltonian = simulator().hamiltonian();
         double dH = newHamiltonian - oldHamiltonian;

         // Compute force bias used in acceptance criterion
         double bias = 0.0;
         for (j = 0; j < nMonomer - 1; ++j) {
            RField<D,T> const & di = dc_[j];
            RField<D,T> const & df = simulator().dc(j);
            computeForceBias(biasField_, di, df, dwc_[j], mobility_);
            bias += Reduce::sum(biasField_);
         }
         bias *= vNode;
         McMove<D,T>::hamiltonianTimer_.stop();

         // Accept or reject move
         McMove<D,T>::decisionTimer_.start();
         bool accept = false;
         double weight = exp(bias - dH);
         accept = McMove<D,T>::random().metropolis(weight);
         if (accept) {
            McMove<D,T>::incrementNAccept();
            simulator().clearState();
         } else {
            simulator().restoreState();
         }
         McMove<D,T>::decisionTimer_.stop();
      }
      McMove<D,T>::totalTimer_.stop();

      return isConverged;
   }

   /*
   * Output data to log file (do-nothing implementation).
   */
   template <int D, class T>
   void ForceBiasMoveBase<D,T>::output()
   {}

}
}
#endif
