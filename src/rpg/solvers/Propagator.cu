/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.h"
#include "Block.h"
#include <prdc/field/cuda/RField.h>
#include <pscf/cuda/VecOp.h>
#include <pscf/cuda/Reduce.h>
#include <pscf/mesh/Mesh.h>

#include <rp/solvers/PropagatorBase.tpp>  // base class implementation

namespace Pscf {
namespace Rp {

   using namespace Util;

   // Public member functions

   /*
   * Constructor.
   */
   template <int D>
   Propagator<D, CudaTp<D> >::Propagator()
    : RpPropagatorT()
   {}

   /*
   * Destructor.
   */
   template <int D>
   Propagator<D, CudaTp<D> >::~Propagator()
   {
      dissociateQFields();

      /*
      * The above function dissociates elements of qFields_ from memory
      * owned by qFieldsAll_. Because this destructor body is called 
      * before the destructors for members, disassociation will occur
      * before either qFieldsAll_ or qFields_ is destroyed.  
      */

   }

   /*
   * Allocate memory used by this propagator.
   */
   template <int D>
   void Propagator<D, CudaTp<D> >::allocate(int ns, const Mesh<D>& mesh)
   {
      RpPropagatorT::allocate(ns, mesh);

      const int meshSize = mesh.size();
      IntVec<D> const & meshDimensions = mesh.dimensions();

      // Allocate memory in qFieldsAll_ using value of ns
      qFieldsAll_.allocate(ns * meshSize);

      // Set up array of associated RField<D> arrays
      qFields_.allocate(ns);
      for (int i = 0; i < ns; ++i) {
         qFields_[i].associate(qFieldsAll_, i*meshSize, meshDimensions);
      }
      isAllocated_ = true;

      PropagatorTmplT::setIsSolved(false);
   }

   /*
   * Reallocate memory used by this propagator using new ns value.
   */
   template <int D>
   void Propagator<D, CudaTp<D> >::reallocate(int ns)
   {
      RpPropagatorT::reallocate(ns);

      // Deallocate memory previously used by this propagator.
      dissociateQFields();      // dissociate, nulliify data pointers
      qFields_.deallocate();    // destroy RField<D> objects for slices
      qFieldsAll_.deallocate(); // destroy contiguous data block

      // Store mesh properties
      Mesh<D> const & mesh = RpPropagatorT::mesh();
      const int meshSize = mesh.size();
      IntVec<D> const & meshDimensions = mesh.dimensions();

      // Allocate memory in qFieldsAll_ using new value of ns
      qFieldsAll_.allocate(ns * meshSize);

      // Recreate associations between qFields_ and qFieldsAll_
      qFields_.allocate(ns);
      for (int i = 0; i < ns; ++i) {
         qFields_[i].associate(qFieldsAll_, i*meshSize, meshDimensions);
      }

      PropagatorTmplT::setIsSolved(false);
   }

   // Private member function

   /*
   * Dissociate qFields_ from associated memory blocks in qFieldsAll_.
   *
   * These associations must be destroyed before qFieldsAll_ is
   * de-allocated or destroyed, because it is an error to deallocate
   * a DeviceArray<T> container that is still referred to by one or
   * more other such containers. Associations are kept track of via
   * a private ReferenceCounter owned by qFieldsAll_ (which stores
   * the number of remaining references), and via a data pointer and
   * CountedReference object owned by each element of qFields_ that
   * has a pointer to the associated ReferenceCounter. For each
   * element of qFields_, invoking the dissociate() member function
   * nullifies the array data pointer, sets the array capacity to
   * zero, decrements the number of references in the associated
   * ReferenceCounter owned by qFieldsAll_, and nullifies the pointer
   * to this ReferenceCounter.
   */
   template <int D>
   void Propagator<D, CudaTp<D> >::dissociateQFields()
   {
      if (qFields_.isAllocated()) {
         int ns = qFields_.capacity();
         for (int i = 0; i < ns; ++i) {
            if (qFields_[i].isAssociated()) {
               qFields_[i].dissociate();
            }
         }
      }
   }

}
}

// Explicit Instantiation definitions
namespace Pscf { 
   namespace Rp {
      template class PropagatorBase<1, CudaTp<1> >;
      template class PropagatorBase<2, CudaTp<2> >;
      template class PropagatorBase<3, CudaTp<3> >;
      template class Propagator<1, CudaTp<1> >;
      template class Propagator<2, CudaTp<2> >;
      template class Propagator<3, CudaTp<3> >;
   }
}
