#ifndef RPG_REAL_MOVE_TPP
#define RPG_REAL_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RealMove.h"
#include "McMove.h" 
#include <rpg/fts/montecarlo/McSimulator.h>
#include <rpg/fts/VecOpFts.h>
#include <prdc/cuda/VecOp.h>
#include <pscf/math/IntVec.h>
#include <util/param/ParamComposite.h>
#include <rpg/System.h>
#include <pscf/cuda/CudaRandom.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   RealMove<D>::RealMove(McSimulator<D>& simulator) 
    : McMove<D>(simulator),
      w_(),
      dwc_(),
      sigma_(0.0),
      isAllocated_(false)
   { setClassName("RealMove"); }

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   RealMove<D>::~RealMove()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void RealMove<D>::readParameters(std::istream &in)
   {
      // Read the probability
      readProbability(in);
      
      // The standard deviation of the Gaussian distribution
      read(in, "sigma", sigma_);
   
   }
   

   template <int D>
   void RealMove<D>::setup()
   {
      McMove<D>::setup();
      const int nMonomer = system().mixture().nMonomer();
      IntVec<D> meshDimensions = system().domain().mesh().dimensions();

      if (!isAllocated_){
         w_.allocate(nMonomer);
         for (int i=0; i < nMonomer; ++i) {
            w_[i].allocate(meshDimensions);
         }
         dwc_.allocate(meshDimensions);
         gaussianField_.allocate(meshDimensions);
         isAllocated_ = true;
      }
   }
   
   /*
   * Attempt unconstrained move
   */
   template <int D>
   void RealMove<D>::attemptMove()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      int i, j;
      
      // Copy current W fields from parent system into w_
      for (i = 0; i < nMonomer; ++i) {
         VecOp::eqV(w_[i], system().w().rgrid(i));
      }
      
      double evec;
      double mean = 0.0;
      
      // Loop over composition eigenvectors of projected chi matrix
      for (j = 0; j < nMonomer - 1; ++j) {
         
         // Generate random numbers between 0.0 and sigma.
         cudaRandom().normal(gaussianField_, sigma_, mean);
         VecOp::eqV(dwc_, gaussianField_);
         
         // Loop over monomer types
         for (i = 0; i < nMonomer; ++i) {
            RField<D> & w = w_[i];
            evec = simulator().chiEvecs(j,i);
            VecOp::addEqVc(w, dwc_, evec);
         }

      }

      // set system r grid
      system().setWRGrid(w_);

   }


   /*
   * Trivial default implementation - do nothing
   */
   template <int D>
   void RealMove<D>::output()
   {}
   
   template<int D>
   void RealMove<D>::outputTimers(std::ostream& out)
   {
      out << "\n";
      out << "RealMove time contributions:\n";
      McMove<D>::outputTimers(out);
   }

}
}
#endif
