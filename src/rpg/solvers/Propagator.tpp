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
         qFields_[i].associate(qFieldsAll_, i*meshSize, meshPtr_->dimensions());
      }

      setIsSolved(false);
   }

   /*
   * Compute initial head QField from final tail QFields of sources.
   */
   template <int D>
   void Propagator<D>::computeHeadThread()
   {

      // Initialize head field (s=0) to 1.0 at all grid points
      int nx = meshPtr_->size();
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
   * Compute initial head QField from final tail QFields of sources.
   */
   template <int D>
   void Propagator<D>::computeHeadBead()
   {
      UTIL_CHECK(PolymerModel::isBead());

      // Set head slice to product to source tail slices
      computeHeadThread();

      // If propagator owns the head vertex, apply the bead field weight
      if (ownsHead()) {
         QField& qh = qFields_[0]; // Head slice of this propagator
         block().stepFieldBead(qh);
      }

   }

   /*
   * Solve the modified diffusion equation for this block.
   */
   template <int D>
   void Propagator<D>::solve()
   {
      UTIL_CHECK(isAllocated());
      
      if (PolymerModel::isThread()) {

         // Initialize head as pointwise product of source propagators
         computeHeadThread();

         // MDE step loop for thread model
         for (int iStep = 0; iStep < ns_ - 1; ++iStep) {
            block().stepThread(qFields_[iStep], qFields_[iStep + 1]);
         }

      } else
      if (PolymerModel::isBead()) {

         // Initialize head, starting with product of source propagators
         computeHeadBead();

         // MDE step loop for bead model (stop before tail vertex)
         int iStep;
         for (iStep = 0; iStep < ns_ - 2; ++iStep) {
            block().stepBead(qFields_[iStep], qFields_[iStep + 1]);
         }

         // Compute q for the tail vertex
         iStep = ns_ - 2;
         if (ownsTail()) {
            // Full step, including bead weight for tail vertex
            block().stepBead(qFields_[iStep], qFields_[iStep + 1]);
         } else {
            // Bond operator, excluding bead weight for tail vertex
            block().stepBondBead(qFields_[iStep], qFields_[iStep + 1]);
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

      // Initialize initial (head) field
      VecOp::eqV(qFields_[0], head);

      for (int iStep = 0; iStep < ns_ - 1; ++iStep) {
         block().stepThread(qFields_[iStep], qFields_[iStep+1]);
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
      if (!isSolved()) {
         UTIL_THROW("Propagator is not solved.");
      }
      if (!hasPartner()) {
         UTIL_THROW("Propagator has no partner set.");
      }
      if (!partner().isSolved()) {
         UTIL_THROW("Partner propagator is not solved");
      }

      QField const & qh = head();
      QField const & qt = partner().tail();
      int nx = meshPtr_->size();
      UTIL_CHECK(nx == qh.capacity());
      UTIL_CHECK(nx == qt.capacity());

      // Compute average product of head slice and partner tail slice
      double Q = 1.0;
      if (PolymerModel::isThread()) {
         Q = block().averageProduct(qh, qt);
      } else {
         UTIL_CHECK(ownsHead() == partner().ownsTail());
         if (ownsHead()) {
            Q = block().averageProductBead(qh, qt);
         } else {
            Q = block().averageProduct(qh, qt);
         }
      }

      return Q;
   }

}
}
#endif
