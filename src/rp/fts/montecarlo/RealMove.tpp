#ifndef RP_REAL_MOVE_TPP
#define RP_REAL_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RealMove.h"
#include <rp/fts/montecarlo/McSimulator.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <pscf/math/IntVec.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   RealMove<D,T>::RealMove(McSimulator<D,T>& simulator)
    : McMoveT(simulator),
      w_(),
      dwc_(),
      sigma_(0.0),
      isAllocated_(false)
   {  ParamComposite::setClassName("RealMove"); }

   /*
   * Read body of parameter file block.
   */
   template <int D, class T>
   void RealMove<D,T>::readParameters(std::istream &in)
   {
      McMoveT::readProbability(in);

      // Standard deviation of field changes
      ParamComposite::read(in, "sigma", sigma_);
   }

   /*
   * Setup before simulation loop.
   */
   template <int D, class T>
   void RealMove<D,T>::setup()
   {
      McMoveT::setup();

      if (!isAllocated_){
         const int nMonomer = system().mixture().nMonomer();
         IntVec<D> meshDimensions = system().domain().mesh().dimensions();
         w_.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            w_[i].allocate(meshDimensions);
         }
         dwc_.allocate(meshDimensions);
         isAllocated_ = true;
      }
   }

   /*
   * Attempt unconstrained move.
   */
   template <int D, class T>
   void RealMove<D,T>::attemptMove()
   {
      // Copy current W fields to w_
      const int nMonomer = system().mixture().nMonomer();
      for (int i = 0; i < nMonomer; ++i) {
         VecOp::eqV(w_[i], system().w().rgrid(i));
      }

      // Loop over composition eigenvectors of projected chi matrix
      double evec, mean;
      mean = 0.0;
      for (int j = 0; j < nMonomer - 1; ++j) {

         // Generate random field changes
         vecRandom().normal(dwc_, sigma_, mean);

         // Add changes to w_ field components
         for (int i = 0; i < nMonomer; ++i) {
            evec = simulator().chiEvecs(j, i);
            VecOp::addEqVc(w_[i], dwc_, evec);
         }
      }

      // Update w-fields of parent system
      system().w().setRGrid(w_);
   }

}
}
#endif
