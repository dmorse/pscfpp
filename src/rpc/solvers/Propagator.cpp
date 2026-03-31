/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.h"
#include "Block.h"
#include <prdc/cpu/RField.h>
#include <pscf/cpu/VecOp.h>
#include <pscf/cpu/Reduce.h>
#include <pscf/mesh/Mesh.h>

#include <rp/solvers/Propagator.tpp>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   Propagator<D>::Propagator()
    : RpPropagatorT()
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
   void Propagator<D>::reallocate(int ns)
   {
      RpPropagatorT::reallocate(ns);
      UTIL_CHECK(RpPropagatorT::ns() == ns);

      // Deallocate all memory previously used by this propagator.
      qFields_.deallocate();

      // NOTE: Variable qFields_ is a DArray< RField<D> > container.
      // The DArray::deallocate() function calls "delete [] ptr", where 
      // ptr is a pointer to the underlying C array. The C++ delete [] 
      // command calls the destructor for each RField<D> array element
      // before deleting the parent array. The RField<D> destructor 
      // deletes the double* array that stores the field associated 
      // with each slice of the propagator. All memory is thus released.

      // Allocate new memory for qFields_ using the new value of ns
      qFields_.allocate(ns);
      for (int i = 0; i < ns; ++i) {
         qFields_[i].allocate(RpPropagatorT::mesh().dimensions());
      }

      PropagatorTmplT::setIsSolved(false);
   }

} // namespace Rpc
} // namespace Pscf

// Explicit instantiation definitions
namespace Pscf { 
   namespace Rp {
      template class Propagator<1, Rpc::Types<1> >;
      template class Propagator<2, Rpc::Types<2> >;
      template class Propagator<3, Rpc::Types<3> >;
   }
   namespace Rpc {
      template class Propagator<1>;
      template class Propagator<2>;
      template class Propagator<3>;
   }
}
