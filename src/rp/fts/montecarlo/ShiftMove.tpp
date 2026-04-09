#ifndef RP_SHIFT_TPP
#define RP_SHIFT_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ShiftMove.h"

#include <pscf/mesh/Mesh.h>
#include <util/containers/Array.h>
#include <util/random/Random.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   ShiftMove<D,T>::ShiftMove(typename T::McSimulator& simulator)
    : McMoveT(simulator),
      maxShift_(0),
      isAllocated_(false)
   {  ParamComposite::setClassName("ShiftMove"); }

   /*
   * Read body of parameter file block.
   */
   template <int D, class T>
   void ShiftMove<D,T>::readParameters(std::istream &in)
   {

      // Read the probability
      McMoveT::readProbability(in);

      // Read the maximum shift
      ParamComposite::read(in, "maxShift", maxShift_);

      // Validate maxShift_ value
      UTIL_CHECK(maxShift_ > 0);
      IntVec<D> const & dim = system().domain().mesh().dimensions();
      for (int i = 0; i < D; i++) {
         UTIL_CHECK(maxShift_ < dim[i]);
      }
   }

   /*
   * Setup just before beginning a simulation.
   */
   template <int D, class T>
   void ShiftMove<D,T>::setup()
   {
      // Setup base class
      McMoveT::setup();

      // Allocate memory if necessary
      if (!isAllocated_) {
         const int nMonomer = system().mixture().nMonomer();
         IntVec<D> const & dim = system().domain().mesh().dimensions();
         w_.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            w_[i].allocate(dim);
         }
         isAllocated_ = true;
      }
   }

   /*
   * Attempt unconstrained move
   */
   template <int D, class T>
   void ShiftMove<D,T>::attemptMove()
   {
      // Select random displacement by integer numbers of mesh points
      IntVec<D> shift;
      for (int i = 0; i < D; i++){
         shift[i] = McMoveT::random().uniformInt(-maxShift_, maxShift_ + 1);
      }

      // Compute shifted fields stored in w_ array
      shiftFields(shift);

      // Update w-fields in parent system
      system().w().setRGrid(w_);
   }

   /*
   * Shift a single field.
   */
   template <int D, class T>
   void ShiftMove<D,T>::shiftField(Array<double> & out, 
                                 Array<double> const & in,
                                 IntVec<D> shift, 
                                 IntVec<D> dimensions) const
   {
      Mesh<D> mesh(dimensions);
      IntVec<D> inPosition, outPosition;
      int meshSize = mesh.size();
      UTIL_CHECK(out.capacity() == meshSize);
      UTIL_CHECK(in.capacity() == meshSize);

      for (int i = 0; i < meshSize; ++i) {
         inPosition = mesh.position(i);
         for (int d = 0; d < D; ++d){
            outPosition[d] = inPosition[d] + shift[d];
            if (outPosition[d] < 0){
               outPosition[d] += dimensions[d];
            }
            UTIL_CHECK(outPosition[d] >= 0);
            outPosition[d] = outPosition[d] % dimensions[d];
         }
         out[mesh.rank(outPosition)] = in[mesh.rank(inPosition)];
      }
   }

}
}
#endif
