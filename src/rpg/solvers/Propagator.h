#ifndef RPG_PROPAGATOR_H
#define RPG_PROPAGATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/Propagator.h>   // base class template
#include <rpg/system/Types.h>        // base class template argument
#include <pscf/cuda/DeviceArray.h>   // member

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;

   /**
   * MDE solver for one direction of one block.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of base class template Rp::Propagator, 
   * and inherit most of their public interface from this base class.  
   * Only a few functions that involve memory allocation on the GPU are 
   * defined or re-defined in this template
   *
   * \see Rp::Propagator
   * \see Pscf::PropagatorTmpl
   * \ingroup Rpg_Solver_Module
   */
   template <int D>
   class Propagator : public Rp::Propagator< D, Types<D> >
   {

   public:

      // Member functions

      /**
      * Constructor.
      */
      Propagator();

      /**
      * Destructor.
      */
      ~Propagator();

      /**
      * Allocate memory used by this propagator.
      *
      * The parameter ns is the number of values of s at which q(r,s) is
      * calculated, including the end values at the terminating vertices.
      * See docs for the function ns(), which returns this value.
      *
      * The address of the associated Mesh<D> object is retained.
      *
      * An Exception is thrown if the propagator is already allocated.
      *
      * \param ns  number of slices (including end points at vertices)
      * \param mesh  spatial discretization mesh
      */
      void allocate(int ns, const Mesh<D>& mesh) override;

      /**
      * Reallocate memory used by this propagator.
      *
      * This function is used when the value of ns is changed, which can
      * occur during some parameter sweeps. See docs for allocate and ns.
      *
      * An Exception is thrown if the propagator has not been previously
      * allocated, or if it is allocated but the value of ns is unchanged.
      *
      * \param ns  number of slices (including end points)
      */
      void reallocate(int ns) override;

      /**
      * Return the full array of q-fields as an unrolled 1D array.
      */
      DeviceArray<cudaReal> const & qAll();

   protected:

      /// Direct base class.
      using RpPropagatorT = Rp::Propagator<D, Types<D> >;

      // Inherited typename alias
      using typename RpPropagatorT::PropagatorTmplT;

      using RpPropagatorT::qFields_;
      using RpPropagatorT::ns_;
      using RpPropagatorT::isAllocated_;
      using RpPropagatorT::computeHead;

   private:

      /**
      * Array containing the entire propagator, stored on the device.
      *
      * The propagator data is stored contiguously to allow batched FFTs
      * to be performed on all contour steps simultaneously, which occurs
      * in Block::computeStress.
      * 
      * Each element of the qFields_ container is an RField<D> that acts
      * as a reference array that points to a slice of the contiguous
      * array qFieldsAll_.  This association is created in the allocate
      * or de-allocate functions, and destroyed by the dissociateQFields
      * function.
      */
      DeviceArray<cudaReal> qFieldsAll_;

      /**
      * Dissociate all qFields_ from associated slices of qFieldsAll_.
      */
      void dissociateQFields();

   };

   // Inline function

   /*
   * Return the full array of q-fields.
   */
   template <int D> inline
   DeviceArray<cudaReal> const & Propagator<D>::qAll()
   {
      UTIL_CHECK(PropagatorTmplT::isSolved());
      return qFieldsAll_;
   }

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class Propagator<1, Rpg::Types<1> >;
      extern template class Propagator<2, Rpg::Types<2> >;
      extern template class Propagator<3, Rpg::Types<3> >;
   }
   namespace Rpg {
      extern template class Propagator<1>;
      extern template class Propagator<2>;
      extern template class Propagator<3>;
   }
}
#endif
