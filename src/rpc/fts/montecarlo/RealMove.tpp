#ifndef RPC_REAL_MOVE_TPP
#define RPC_REAL_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RealMove.h"
#include "McMove.h" 
#include <rpc/fts/montecarlo/McSimulator.h>
#include <util/param/ParamComposite.h>
#include <rpc/System.h>
#include <util/random/Random.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   RealMove<D>::RealMove(McSimulator<D>& simulator) 
    : McMove<D>(simulator),
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
      IntVec<D> const & meshDimensions = system().domain().mesh().dimensions();
      if (!isAllocated_){
         dwc_.allocate(meshDimensions);
         w_.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            w_[i].allocate(meshDimensions);
         }
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
      
      // Copy current fields to w_
      for (int i = 0; i < nMonomer; ++i) {
         w_[i] = system().w().rgrid(i);
      }
      
      // Loop over composition eigenvectors of projected chi matrix
      for (int j = 0; j < nMonomer - 1; j++){

         // Generate Gaussian distributed random numbers
         for (int k = 0; k < meshSize; k++){
            dwc_[k] = sigma_* random().gaussian();
         }
         
         // Loop over monomer types
         double evec;
         for (int i = 0; i < nMonomer; ++i) {
            RField<D> & w = w_[i];
            evec = simulator().chiEvecs(j, i);
            for (int k = 0; k < meshSize; ++k) {
               w[k] += evec*dwc_[k];
            }
         }

      }

      // Update w-fields in parent system
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
