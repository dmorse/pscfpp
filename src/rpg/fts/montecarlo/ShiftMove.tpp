#ifndef RPG_SHIFT_TPP
#define RPG_SHIFT_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ShiftMove.h"
#include "McMove.h"
#include <rpg/solvers/Mixture.h>
#include <rpg/field/Domain.h>
#include <rpg/fts/montecarlo/McSimulator.h>
#include <rpg/system/System.h>
#include <rpg/fts/VecOpFts.h>
#include <pscf/cuda/VecOp.h>
#include <pscf/math/IntVec.h>
#include <pscf/cuda/CudaVecRandom.h>
#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <util/param/ParamComposite.h>
#include <util/random/Random.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   ShiftMove<D>::ShiftMove(McSimulator<D>& simulator)
    : McMove<D>(simulator),
      maxShift_(0),
      isAllocated_(false)
   { setClassName("ShiftMove"); }

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   ShiftMove<D>::~ShiftMove()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void ShiftMove<D>::readParameters(std::istream &in)
   {

      // Read the probability
      readProbability(in);

      // The standard deviation of the Gaussian distribution
      read(in, "maxShift", maxShift_);
      IntVec<D> const & meshDimensions = system().domain().mesh().dimensions();
      for (int i = 0; i < D; i++){
         UTIL_CHECK(maxShift_ < meshDimensions[i]);
      }
   }

   template <int D>
   void ShiftMove<D>::setup()
   {
      McMove<D>::setup();
      const int nMonomer = system().mixture().nMonomer();
      IntVec<D> const & meshDimensions = system().domain().mesh().dimensions();

      if (!isAllocated_){
         w0_.allocate(nMonomer);
         w_.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            w0_[i].allocate(meshDimensions);
            w_[i].allocate(meshDimensions);
         }
         isAllocated_ = true;
      }
   }

   /*
   * Attempt unconstrained move
   */
   template <int D>
   void ShiftMove<D>::attemptMove()
   {
      const int nMonomer = system().mixture().nMonomer();
      IntVec<D> const & meshDimensions = system().domain().mesh().dimensions();

      // Copy current W fields from parent system into w0_
      for (int i = 0; i < nMonomer; ++i) {
         VecOp::eqV(w0_[i], system().w().rgrid(i));
      }

      // Shift Field direction
      IntVec<D> shift;
      IntVec<D> position;
      IntVec<D> shiftPosition;
      Mesh<D> mesh(meshDimensions);
      MeshIterator<D> iter(meshDimensions);
      for (int i = 0; i < D; i++){
         shift[i] = simulator().random().uniformInt(-maxShift_, maxShift_ + 1);
      }

      for (int j = 0; j< nMonomer; ++j){
         VecOpFts::shiftWField(w_[j], w0_[j], meshDimensions, shift);
      }

      // Update w-fields in parent system
      system().w().setRGrid(w_);
   }

   /*
   * Trivial default implementation - do nothing
   */
   template <int D>
   void ShiftMove<D>::output()
   {}

   template<int D>
   void ShiftMove<D>::outputTimers(std::ostream& out)
   {
      out << "\n";
      out << "ShiftMove time contributions:\n";
      McMove<D>::outputTimers(out);
   }

}
}
#endif
