#ifndef RPG_PROPAGATOR_TPP
#define RPG_PROPAGATOR_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.h"
#include "Block.h"
#include <pscf/cuda/GpuResources.h>
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
   void Propagator<D>::computeHead()
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
               UTIL_THROW("Source not solved in computeHead");
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
      UTIL_CHECK(isAllocated());

      computeHead();
      for (int iStep = 0; iStep < ns_ - 1; ++iStep) {
         block().step(qFields_[iStep], qFields_[iStep+1]);
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
         block().step(qFields_[iStep], qFields_[iStep+1]);
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

      // Take inner product of head and partner tail fields
      // cannot reduce assuming one propagator, qh == 1
      // polymers are divided into blocks midway through
      int nx = meshPtr_->size();
      double Q = Reduce::innerProduct(head(), partner().tail()) / double(nx);
      return Q;
   }

}
}
#endif
