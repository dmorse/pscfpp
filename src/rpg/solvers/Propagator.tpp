#ifndef RPG_PROPAGATOR_TPP
#define RPG_PROPAGATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.h"
#include "Block.h"
#include <prdc/cuda/resources.h>
#include <pscf/mesh/Mesh.h>
#include <util/containers/DArray.h>

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   Propagator<D>::Propagator()
    : blockPtr_(0),
      meshPtr_(0),
      ns_(0),
      isAllocated_(false)
   {}

   /*
   * Destructor.
   */
   template <int D>
   Propagator<D>::~Propagator()
   {}

   /*
   * Allocate memory used by this propagator.
   */
   template <int D>
   void Propagator<D>::allocate(int ns, const Mesh<D>& mesh)
   {
      ns_ = ns;
      meshPtr_ = &mesh;

      // Allocate memory in qFieldsAll_ using value of ns
      int meshSize = meshPtr_->size();
      qFieldsAll_.allocate(ns * meshSize);

      // Set up array of associated RField<D> arrays
      qFields_.allocate(ns);
      for (int i = 0; i < ns; ++i) {
         qFields_[i].associate(qFieldsAll_, i*meshSize, meshPtr_->dimensions());
      }
      isAllocated_ = true;
   }

   /*
   * Reallocate memory used by this propagator using new ns value.
   */
   template <int D>
   void Propagator<D>::reallocate(int ns)
   {
      // Preconditions
      UTIL_CHECK(meshPtr_);
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(ns_ != ns);
      ns_ = ns;

      // Deallocate memory previously used by this propagator.
      qFields_.deallocate(); // Destroys associated RField<D> objects 
                             // but not the underlying data
      qFieldsAll_.deallocate(); // Destroy actual propagator data

      // Allocate memory in qFieldsAll_ using new value of ns
      int meshSize = meshPtr_->size();
      qFieldsAll_.allocate(ns * meshSize);

      // Set up array of associated RField<D> arrays
      qFields_.allocate(ns);
      for (int i = 0; i < ns; ++i) {
         qFields_[i].associate(qFieldsAll_, i*meshSize, 
                               meshPtr_->dimensions());
      }

      setIsSolved(false);
   }

   /*
   * Compute initial head q-field from final tail q-fields of sources.
   */
   template <int D>
   void Propagator<D>::computeHead()
   {
      UTIL_CHECK(meshPtr_);
      int nx = meshPtr_->size();

      // Initialize head field (s=0) to 1.0 at all grid points
      VecOp::eqS(qFields_[0], 1.0);

      // Multiply head q-field by tail q-fields of all sources
      if (nSource() > 0) {
         DArray<DeviceArray<cudaReal> const *> tails;
         tails.allocate(nSource()+1);
         tails[0] = &qFields_[0]; 
         for (int is = 0; is < nSource(); ++is) {
            if (!source(is).isSolved()) {
               UTIL_THROW("Source not solved in computeHeadThread");
            }
            tails[is+1] = &(source(is).tail());
         }
         VecOp::mulVMany(qFields_[0], tails);
      }
   }

   /*
   * Solve the modified diffusion equation for this block.
   */
   template <int D>
   void Propagator<D>::solve()
   {
      UTIL_CHECK(blockPtr_);
      UTIL_CHECK(isAllocated());
      
      // Initialize head, starting with product of source propagators
      computeHead();

      if (PolymerModel::isThread()) {

         // MDE step loop for thread model
         for (int iStep = 0; iStep < ns_ - 1; ++iStep) {
            block().stepThread(qFields_[iStep], qFields_[iStep + 1]);
         }

      } else
      if (PolymerModel::isBead()) {

         // Half-bond and bead weight for first bead
         if (isHeadEnd()) {
            VecOp::eqV(qFields_[1], qFields_[0]);
         } else {
            block().stepHalfBondBead(qFields_[0], qFields_[1]);
         }
         block().stepFieldBead(qFields_[1]);

         // MDE step loop for bead model (stop before tail vertex)
         int iStep;
         for (iStep = 1; iStep < ns_ - 2; ++iStep) {
            block().stepBead(qFields_[iStep], qFields_[iStep + 1]);
         }

         // Half-bond for tail slice
         if (isTailEnd()) {
            VecOp::eqV(qFields_[ns_-2], qFields_[ns_-2]);
         } else {
            block().stepHalfBondBead(qFields_[ns_-2], qFields_[ns_-1]);
         }

      } else {
         // This should be impossible
         UTIL_THROW("Unexpected PolymerModel type");
      }

      setIsSolved(true);
   }

   /*
   * Solve the modified diffusion equation with specified initial field.
   */
   template <int D>
   void Propagator<D>::solve(RField<D> const & head)
   {
      UTIL_CHECK(isAllocated());

      // Initialize head slice (index 0)
      VecOp::eqV(qFields_[0], head);

      if (PolymerModel::isThread()) {

         // MDE step loop for thread model
         for (int iStep = 0; iStep < ns_ - 1; ++iStep) {
            block().stepThread(qFields_[iStep], qFields_[iStep + 1]);
         }

      } else
      if (PolymerModel::isBead()) {

         // Half-bond and bead weight for first bead
         if (isHeadEnd()) {
            VecOp::eqV(qFields_[1], qFields_[0]);
         } else {
            block().stepHalfBondBead(qFields_[0], qFields_[1]);
         }
         block().stepFieldBead(qFields_[1]);

         // MDE step loop for bead model (stop before tail vertex)
         int iStep;
         for (iStep = 1; iStep < ns_ - 2; ++iStep) {
            block().stepBead(qFields_[iStep], qFields_[iStep + 1]);
         }

         // Half-bond for tail slice
         if (!isTailEnd()) {
            VecOp::eqV(qFields_[ns_-2], qFields_[ns_-2]);
         } else {
            block().stepHalfBondBead(qFields_[ns_-2], qFields_[ns_-1]);
         }

      } else {
         // This should be impossible
         UTIL_THROW("Unexpected PolymerModel type");
      }
      setIsSolved(true);
   }

   /*
   * Integrate to calculate monomer concentration for this block
   */
   template <int D>
   double Propagator<D>::computeQ()
   {
      // Preconditions
      UTIL_CHECK(meshPtr_);
      UTIL_CHECK(isAllocated_);
      if (!isSolved()) {
         UTIL_THROW("Propagator is not solved.");
      }
      if (!hasPartner()) {
         UTIL_THROW("Propagator has no partner set.");
      }
      if (!partner().isSolved()) {
         UTIL_THROW("Partner propagator is not solved");
      }
      UTIL_CHECK(isHeadEnd() == partner().isTailEnd());

      double Q = 0.0;
      if (PolymerModel::isBead() && isHeadEnd()) {
         // Compute average of q for last bead of partner
         RField<D> const & qt = partner().q(ns_-2);
         Q = Reduce::sum(qt);
      } else {
         // Compute average product of head slice and partner tail slice
         RField<D> const & qh = head();
         RField<D> const & qt = partner().tail();
         Q = Reduce::innerProduct(qh, qt);
      }
      Q /= double(meshPtr_->size());
      return Q;
   }

}
}
#endif

