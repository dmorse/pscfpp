#ifndef PSPC_PROPAGATOR_TPP
#define PSPC_PROPAGATOR_TPP

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
namespace Pspc {

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
   }

   /*
   * Compute initial head QField from final tail QFields of sources.
   */
   template <int D>
   void Propagator<D>::computeHead()
   {

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
   * Solve the modified diffusion equation for this block.
   */
   template <int D>
   void Propagator<D>::solve()
   {
      UTIL_CHECK(isAllocated());
      computeHead();
      for (int iStep = 0; iStep < ns_ - 1; ++iStep) {
         block().step(qFields_[iStep], qFields_[iStep + 1]);
      }
      setIsSolved(true);
   }

   /*
   * Solve the modified diffusion equation with a specified initial condition.
   */
   template <int D>
   void Propagator<D>::solve(QField const & head)
   {
      int nx = meshPtr_->size();
      UTIL_CHECK(head.capacity() == nx);

      // Initialize initial (head) field
      QField& qh = qFields_[0];
      for (int i = 0; i < nx; ++i) {
         qh[i] = head[i];
      }

      // Setup solver and solve
      for (int iStep = 0; iStep < ns_ - 1; ++iStep) {
         block().step(qFields_[iStep], qFields_[iStep + 1]);
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
      QField const& qh = head();
      QField const& qt = partner().tail();
      int nx = meshPtr_->size();
      UTIL_CHECK(qt.capacity() == nx);
      UTIL_CHECK(qh.capacity() == nx);

      // Take inner product of head and partner tail fields
      double Q = 0;
      for (int i =0; i < nx; ++i) {
         Q += qh[i]*qt[i];
      }
      Q /= double(nx);
      return Q;
   }

}
}
#endif
