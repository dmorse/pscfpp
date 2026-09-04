#ifndef RP_PROPAGATOR_C_H
#define RP_PROPAGATOR_C_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/solvers/PropagatorBase.h>  // base class template
#include <pscf/backend/CPT.h>          // specialized argument

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Pscf::Prdc;

   // Declare primary template
   template <int D, class T> class Propagator;

   /**
   * MDE solver for one direction of one block (CPU specialization).
   *
   * This is the partial specialization for use with a serial CPU backend. 
   * This specialization inherits its public interface and almost all of
   * its source code from the PropagatorBase base class template. 
   *
   * Virtual allocate() and reallocate() functions are defined here rather
   * than in PropagatorBase because the CPU and  GPU versions of this 
   * class use different strategies for memory allocation.
   *
   * \see PropagatorBase
   * \see Pscf::PropagatorTmpl
   * \ingroup Rp_Solver_Module
   */
   template <int D>
   class Propagator<D,CPT> 
    : public PropagatorBase<D,CPT>
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

   protected:

      /// Direct (parent) base class.
      using RpPropagatorT = PropagatorBase<D,CPT>;

      /// Inherited typename alias for indirect (grandparent) base class.
      using typename RpPropagatorT::PropagatorTmplT;

      /// Inherited non-dependent class members.
      using RpPropagatorT::qFields_;
      using RpPropagatorT::ns_;
      using RpPropagatorT::isAllocated_;
      using RpPropagatorT::computeHead;

   };

} // namespace Rp
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class PropagatorBase<1,CPT>;
      extern template class PropagatorBase<2,CPT>;
      extern template class PropagatorBase<3,CPT>;
      extern template class Propagator<1,CPT>;
      extern template class Propagator<2,CPT>;
      extern template class Propagator<3,CPT>;
   }
}
#endif
