#ifndef RPC_PROPAGATOR_TPP
#define RPC_PROPAGATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.h"
#include "Block.h"

#include <pscf/mesh/Mesh.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   Propagator<D>::Propagator()
    : blockPtr_(nullptr),
      meshPtr_(nullptr),
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
      UTIL_CHECK(!isAllocated_);
      ns_ = ns;
      meshPtr_ = &mesh;

      qFields_.allocate(ns);
      for (int i = 0; i < ns; ++i) {
         qFields_[i].allocate(mesh.dimensions());
      }
      isAllocated_ = true;
   }

   /*
   * Reallocate memory used by this propagator using new ns value.
   */
   template <int D>
   void Propagator<D>::reallocate(int ns)
   {
      UTIL_CHECK(meshPtr_);
      UTIL_CHECK(isAllocated_);
      UTIL_CHECK(ns_ != ns);
      ns_ = ns;

      // Deallocate all memory previously used by this propagator.
      qFields_.deallocate();

      // NOTE: The qFields_ container is a DArray<QField>, where QField
      // is a typedef for DFields<D>. The DArray::deallocate() function
      // calls "delete [] ptr", where ptr is a pointer to the underlying
      // C array. The C++ delete [] command calls the destructor for
      // each element in an array before deleting the array itself.
      // The RField<D> destructor deletes the the double* array that 
      // stores the field associated with each slice of the propagator.

      // Allocate new memory for qFields_ using new value of ns
      qFields_.allocate(ns);
      for (int i = 0; i < ns; ++i) {
         qFields_[i].allocate(meshPtr_->dimensions());
      }

      setIsSolved(false);
   }

   /*
   * Compute initial head QField for the thread model.
   */
   template <int D>
   void Propagator<D>::computeHeadThread()
   {
      UTIL_CHECK(meshPtr_);

      // Reference to head of this propagator
      QField& qh = qFields_[0];

      // Initialize qh field to 1.0 at all grid points
      int ix;
      int nx = meshPtr_->size();
      for (ix = 0; ix < nx; ++ix) {
         qh[ix] = 1.0;
      }

      // Pointwise multiply tail QFields of all sources
      for (int is = 0; is < nSource(); ++is) {
         if (!source(is).isSolved()) {
            UTIL_THROW("Source not solved in computeHead");
         }
         QField const& qt = source(is).tail();
         for (ix = 0; ix < nx; ++ix) {
            qh[ix] *= qt[ix];
         }
      }
   }

   /*
   * Compute initial head QField for the bead model.
   */
   template <int D>
   void Propagator<D>::computeHeadBead()
   {
      UTIL_CHECK(blockPtr_);
      UTIL_CHECK(PolymerModel::isBead());

      // Set head slice to product of source tail slices
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
      UTIL_CHECK(blockPtr_);
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

         // Initialize head slice, starting with product of sources
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
   * Solve the MDE with a specified initial condition at the head.
   */
   template <int D>
   void Propagator<D>::solve(QField const & head)
   {
      UTIL_CHECK(blockPtr_);
      UTIL_CHECK(meshPtr_);
      int nx = meshPtr_->size();
      UTIL_CHECK(head.capacity() == nx);

      // Initialize initial (head) field
      QField& qh = qFields_[0];
      for (int i = 0; i < nx; ++i) {
         qh[i] = head[i];
      }

      // Loop over steps to solve
      if (PolymerModel::isThread()) {

         // MDE step loop for thread model
         for (int iStep = 0; iStep < ns_ - 1; ++iStep) {
            block().stepThread(qFields_[iStep], qFields_[iStep + 1]);
         }

      } else 
      if (PolymerModel::isBead()) {

         // MDE step loop for bead model (step before the tail vertex)
         int iStep;
         for (iStep = 0; iStep < ns_ - 2; ++iStep) {
            block().stepBead(qFields_[iStep], qFields_[iStep + 1]);
         }

         // Compute q for the tail vertex
         iStep = ns_ - 2;
         if (ownsTail()) {
            // Include bead weight exp(W[i]*ds) for tail vertex
            block().stepBead(qFields_[iStep], qFields_[iStep + 1]);
         } else {
            // Exclude bead weight for tail vertex
            block().stepBondBead(qFields_[iStep], qFields_[iStep + 1]);
         }

      } else {
         // This should be impossible
         UTIL_THROW("Unexpected PolymerModel type");
      }

      setIsSolved(true);
   }

   /*
   * Compute spatial average of product of head and tail of partner.
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
      UTIL_CHECK(blockPtr_);
      UTIL_CHECK(meshPtr_);

      QField const& qh = head();
      QField const& qt = partner().tail();
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
