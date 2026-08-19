/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>

#include <pscf/backends/CPT.h>
#include <rp/solvers/PropagatorBase.tpp>
#include "Propagator_c.h"

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   Propagator<D,CPT>::Propagator()
    : RpPropagatorT()
   {}

   /*
   * Destructor.
   */
   template <int D>
   Propagator<D,CPT>::~Propagator()
   {}

   /*
   * Allocate memory used by this propagator.
   */
   template <int D>
   void Propagator<D,CPT>::allocate(int ns, const Mesh<D>& mesh)
   {
      RpPropagatorT::allocate(ns, mesh);
      UTIL_CHECK(RpPropagatorT::ns() == ns);

      qFields_.allocate(ns);
      for (int i = 0; i < ns; ++i) {
         qFields_[i].allocate(mesh.dimensions());
      }
      isAllocated_ = true;

      PropagatorTmplT::setIsSolved(false);
   }

   /*
   * Reallocate memory used by this propagator using new ns value.
   */
   template <int D>
   void Propagator<D,CPT>::reallocate(int ns)
   {
      RpPropagatorT::reallocate(ns);
      UTIL_CHECK(RpPropagatorT::ns() == ns);

      // Deallocate all memory previously used by this propagator.
      qFields_.deallocate();

      // NOTE: Variable qFields_ is a DArray< RField<D,CPT> > container.
      // The DArray::deallocate() function calls "delete [] ptr", where 
      // ptr is a pointer to the underlying C array. The C++ delete [] 
      // command calls the destructor for each RField<D,CPT> array element
      // before deleting the parent array. The RField<D,CPT> destructor 
      // deletes the double* array that stores the field associated 
      // with each slice of the propagator. All memory is thus released.

      // Allocate new memory for qFields_ using the new value of ns
      qFields_.allocate(ns);
      for (int i = 0; i < ns; ++i) {
         qFields_[i].allocate(RpPropagatorT::mesh().dimensions());
      }

      PropagatorTmplT::setIsSolved(false);
   }

} // namespace Rp
} // namespace Pscf

// Explicit instantiation definitions
namespace Pscf { 
   namespace Rp {
      template class PropagatorBase<1,CPT>;
      template class PropagatorBase<2,CPT>;
      template class PropagatorBase<3,CPT>;
      template class Propagator<1,CPT>;
      template class Propagator<2,CPT>;
      template class Propagator<3,CPT>;
   }
}
