#ifndef RPC_SHIFT_TPP
#define RPC_SHIFT_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ShiftMove.h"
#include "McMove.h"
#include <rpc/fts/montecarlo/McSimulator.h>
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <pscf/mesh/Mesh.h>
#include <util/param/ParamComposite.h>
#include <util/random/Random.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   ShiftMove<D>::ShiftMove(Rp::McSimulator<D, Rpc::Types<D> >& simulator)
    : Rp::McMove<D, Rpc::Types<D> >(simulator),
      maxShift_(0),
      isAllocated_(false)
   {  ParamComposite::setClassName("ShiftMove"); }

   /*
   * Read body of parameter file block.
   */
   template <int D>
   void ShiftMove<D>::readParameters(std::istream &in)
   {

      // Read the probability
      Rp::McMove<D, Rpc::Types<D> >::readProbability(in);

      // Read the maximum shift
      ParamComposite::read(in, "maxShift", maxShift_);

      // Validate maxShift_ value
      UTIL_CHECK(maxShift_ > 0);
      IntVec<D> const & dim = system().domain().mesh().dimensions();
      for (int i = 0; i < D; i++){
         UTIL_CHECK(maxShift_ < dim[i]);
      }
   }

   template <int D>
   void ShiftMove<D>::setup()
   {
      // Setup base class
      Rp::McMove<D, Rpc::Types<D> >::setup();

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
   template <int D>
   void ShiftMove<D>::attemptMove()
   {
      // Select random shift
      IntVec<D> shift;
      for (int i = 0; i < D; i++){
         shift[i] = random().uniformInt(-maxShift_, maxShift_ + 1);
      }

      // Compute shifted fields stored in w_ array
      shiftFields(shift);

      // Update w-fields in parent system
      system().w().setRGrid(w_);
   }

   /*
   * Compute and store array w_ of shifted fields.
   */
   template <int D>
   void ShiftMove<D>::shiftFields(IntVec<D> const & shift)
   {
      IntVec<D> const & dimensions = system().domain().mesh().dimensions();
      const int nMonomer = system().mixture().nMonomer();
      for (int j = 0; j< nMonomer; ++j) {
         RField<D> const & w0 = system().w().rgrid(j);
         shiftField(w_[j], w0, shift, dimensions);
      }
   }

   /*
   * Shift a single field.
   */
   template <int D>
   void ShiftMove<D>::shiftField(Array<double> & out, 
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
      Rp::McMove<D, Rpc::Types<D> >::outputTimers(out);
   }

}
}
#endif
