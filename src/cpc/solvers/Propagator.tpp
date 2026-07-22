#ifndef CPC_PROPAGATOR_TPP
#define CPC_PROPAGATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.h"
#include "Block.h"

#include <pscf/cpu/complex.h>
#include <pscf/math/arithmetic.h>
#include <pscf/mesh/Mesh.h>

namespace Pscf {
namespace Cpc {

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
   * Associate this propagator with a block.
   */
   template <int D> inline 
   void Propagator<D>::setBlock(Block<D>& block)
   {
      UTIL_CHECK(!blockPtr_);  
      blockPtr_ = &block; 
   }

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
   * Reallocate memory used by this propagator using a new ns value.
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

      // NOTE: The qFields_ container is a DArray<FieldT>, where FieldT
      // is an alias for CField<D,CPT>. The DArray::deallocate() function
      // calls "delete [] ptr", where ptr is a pointer to the underlying
      // C array. The C++ delete [] command calls the destructor for each
      // element of type CField<T> before deleting the outer array itself. 
      // The CField<D,CPT> destructor deletes the fftw_complex* array that 
      // stores the field associated with each slice of the propagator.

      // Allocate new memory for qFields_ using the new value of ns
      qFields_.allocate(ns);
      for (int i = 0; i < ns; ++i) {
         qFields_[i].allocate(meshPtr_->dimensions());
      }

      setIsSolved(false);
   }

   /*
   * Compute initial head q-field for the thread model.
   */
   template <int D>
   void Propagator<D>::computeHead()
   {
      UTIL_CHECK(meshPtr_);
      int nx = meshPtr_->size();

      // Reference to head of this propagator
      FieldT& qh = qFields_[0];

      // Initialize qh field to 1.0 at all grid points
      for (int ix = 0; ix < nx; ++ix) {
         assign(qh[ix], 1.0);
         //qh[ix] = 1.0;
      }
      // VecOp::eqS(qh, 1.0);

      if (!isHeadEnd()) {
         // Pointwise multiply tail q-fields of all sources
         for (int is = 0; is < nSource(); ++is) {
            if (!source(is).isSolved()) {
               UTIL_THROW("Source not solved in computeHead");
            }
            FieldT const& qt = source(is).tail();
            for (int ix = 0; ix < nx; ++ix) {
               mulEq(qh[ix], qt[ix]);
            }
            // VecOp::mulEq(qh, qt);
         }
      }

   }

   /*
   * Assign one field to another.
   */
   template <int D>
   void Propagator<D>::assignField(FieldT& lhs, FieldT const & rhs)
   {
      int nx = lhs.capacity();
      UTIL_CHECK(rhs.capacity() == nx);
      for (int ix = 0; ix < nx; ++ix) {
          assign(lhs[ix], rhs[ix]);
          // lhs[ix] = rhs[ix];
      }
      // VecOp::eqV(lhs, rhs);
   }

   /*
   * Solve the modified diffusion equation for this block.
   */
   template <int D>
   void Propagator<D>::solve()
   {
      UTIL_CHECK(blockPtr_);
      UTIL_CHECK(isAllocated());

      // Initialize head as pointwise product of source propagators
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
            assignField(qFields_[1], qFields_[0]);
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
            assignField(qFields_[ns_-1], qFields_[ns_-2]);
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
   * Solve the MDE with a specified initial condition at the head.
   */
   template <int D>
   void Propagator<D>::solve(FieldT const & head)
   {
      UTIL_CHECK(blockPtr_);
      UTIL_CHECK(meshPtr_);
      UTIL_CHECK(isAllocated_);
      int nx = meshPtr_->size();
      UTIL_CHECK(head.capacity() == nx);

      // Initialize initial (head) field
      FieldT& qh = qFields_[0];
      for (int i = 0; i < nx; ++i) {
         assign(qh[i], head[i]);
      }
      // VecOp::eqV(qh, head);

      if (PolymerModel::isThread()) {

         // MDE step loop for thread model
         for (int iStep = 0; iStep < ns_ - 1; ++iStep) {
            block().stepThread(qFields_[iStep], qFields_[iStep + 1]);
         }

      } else 
      if (PolymerModel::isBead()) {

         // Half-bond and bead weight for first bead
         if (isHeadEnd()) {
            assignField(qFields_[1], qFields_[0]);
         } else {
            block().stepHalfBondBead(qFields_[0], qFields_[1]);
         }
         block().stepFieldBead(qFields_[1]);

         // MDE step loop for bead model (stop before tail vertex)
         int iStep;
         for (iStep = 1; iStep < ns_ - 2; ++iStep) {
            block().stepBead(qFields_[iStep], qFields_[iStep + 1]);
         }

         // Half-bond for tail
         if (isTailEnd()) {
            assignField(qFields_[ns_-1], qFields_[ns_-2]);
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
   * Compute spatial average of product of head and tail of partner.
   */
   template <int D>
   void Propagator<D>::computeQ(std::complex<double> & Q) const
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
      UTIL_CHECK(isHeadEnd() == partner().isTailEnd());
      UTIL_CHECK(meshPtr_);
      int nx = meshPtr_->size();

      fftw_complex sum;
      assign(sum, 0.0);
      if (PolymerModel::isBead() && isHeadEnd()) {
         // Compute average of q for last bead of partner
         FieldT const& qt = partner().q(ns_-2);
         for (int ix = 0; ix < nx; ++ix) {
            addEq(sum, qt[ix]);
         }
      } else {
         // Compute average product of head slice and partner tail slice
         FieldT const& qh = head();
         FieldT const& qt = partner().tail();
         fftw_complex product;
         for (int ix = 0; ix < nx; ++ix) {
            mul(product, qh[ix], qt[ix]);
            addEq(sum, product);
            // sum += qh[ix]*qt[ix];
         }
      }
      double meshSize = double(nx);
      divEq(sum, meshSize);
      assign(Q, sum);
   }

}
}
#endif
