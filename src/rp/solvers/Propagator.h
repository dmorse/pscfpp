#ifndef RP_PROPAGATOR_H
#define RP_PROPAGATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/solvers/PropagatorTmpl.h> // base class template
#include <util/containers/DArray.h>      // member


// Forward declaration
namespace Pscf {
   template <int D> class Mesh;
}

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * MDE solver for one direction of one block.
   *
   * A fully initialized Propagator has associations with a parent Block 
   * object that owns this propagator, and with a partner Propagator that 
   * solves the MDE within the same block in the opposite direction. It 
   * also has an association with a Mesh<D> that describes a spatial grid, 
   * and associations with zero or more source Propagator objects that 
   * are used to compute an initial condition for this propagator at its 
   * head vertex.
   *
   * The parent Block stores information required to numerically
   * solve the modified diffusion equation (MDE), including quantities
   * that depend upon the w-field associated with this block, the unit
   * cell parameters and (in the thread model) the contour step size.
   * These quantities are set and stored by the block because their values
   * are the same for the two propagators owned by each block, but may be
   * different for different blocks.  The algorithm used by a Propagator
   * to solve the MDE repeatedly calls step functions provided by the
   * parent Block.
   *
   * Specializations of this class template are used as base classes for 
   * two closely analogous class templates, also named Propagator, that 
   * are defined in the Rpc and Rpg program-level name spaces for use in 
   * the pscf_rpc and pscf_rpg programs, respectively. 
   *
   * <b> Template parameters </b>:
   *
   *   - D  dimension of space
   *   - T  class with aliases for use in a program-level namespace.
   *
   * \ingroup Rp_Solver_Module
   */
   template <int D, class T>
   class Propagator : public PropagatorTmpl<typename T::Propagator>
   {

   public:

      /**
      * Associate this propagator with a block.
      *
      * \param block associated Block object
      */
      void setBlock(typename T::Block& block);

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
      virtual void allocate(int ns, const Mesh<D>& mesh);

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
      virtual void reallocate(int ns);

      /**
      * Solve the modified diffusion equation (MDE) for this block.
      *
      * This function computes an initial q-field at the head of this
      * block, and then solves the modified diffusion equation (MDE) to
      * propagate the solution from the head to the tail. Algorithms for
      * the thread or bead model may be used, depending on the value of
      * PolymerModel::model().
      */
      void solve();

      /**
      * Solve the MDE for a specified initial condition.
      *
      * This function solves the modified diffusion equation (MDE) for
      * this block with a specified initial condition, which is given by
      * the function parameter "head".
      *
      * \param head  initial condition of q-field at head of block
      */
      void solve(typename T::RField const & head);

      /**
      * Compute and return partition function for the polymer molecule.
      *
      * This function computes the partition function Q for the molecule
      * as a spatial average of the pointwise product of the initial/head
      * slice for this propagator and the final/tail slice of its partner.
      *
      * \param Q  output value, spatial average of q*q^{+} at head
      */
      void computeQ(double & Q) const;

      /**
      * Return q-field at a specified step.
      *
      * \param i  step index, 0 <= i < ns
      */
      const typename T::RField& q(int i) const;

      /**
      * Return q-field at the initial (head) vertex.
      */
      const typename T::RField& head() const;

      /**
      * Return q-field at the terminal (tail) vertex.
      *
      * This function throws an Exception if invoked while the bead model
      * is in use (i.e., if PolymerModel::isThread() == false) and the tail
      * for this propagator is a chain end (i.e., if isTailEnd() == true).
      * In this case, the tail slice is not needed, and so is not computed.
      */
      const typename T::RField& tail() const;

      /**
      * Get the associated Block object by const reference.
      */
      typename T::Block const & block() const;

      /**
      * Return the associated Mesh object by const reference.
      */
      Mesh<D> const & mesh() const;

      /**
      * Get the number of values of s (or slices), including head and tail.
      *
      * The value of ns is the number of values of s at which q(r,s) is
      * calculated, including the end values at the terminating vertices
      * (the head and tail).  In the bead model, this is two more than the
      * number of beads in the block. In the thread model, this is one
      * more than the number length/ds of contour length steps.
      */
      int ns() const;

      /**
      * Has memory been allocated for this propagator?
      */
      bool isAllocated() const;

      /// Direct (parent) base class.
      using PropagatorTmplT = PropagatorTmpl<typename T::Propagator>;

      // Inherited non-dependent members (selected, for convenience)
      using PropagatorTmplT::source;
      using PropagatorTmplT::partner;

   protected:

      /// Array of propagator slices at different contour variable values.
      DArray<typename T::RField> qFields_;

      /// Number of slices, including head and tail slices.
      int ns_;

      /// Is this propagator allocated?
      bool isAllocated_;

      /**
      * Constructor.
      */
      Propagator();

      /**
      * Destructor.
      */
      ~Propagator();

      /**
      * Compute initial q-field at the head vertex.
      *
      * In either the thread or bead model, the head slice of each
      * propagator is the product of tail slices for incoming propagators
      * from other bonds that terminate at the head vertex, or 1 for a
      * chain end.
      */
      void computeHead();

      /**
      * Compute solution of modified diffusion equation (MDE).
      *
      * This function accesses a precomputed head slice qFields_[0]
      * and computes the rest of the propagator.
      */
      void solveMde();

      /**
      * Get the associated Block object by non-const reference.
      */
      typename T::Block& block();

   private:

      /// Pointer to the associated Block.
      typename T::Block* blockPtr_;

      /// Pointer to the associated Mesh.
      Mesh<D> const * meshPtr_;

   };

   // Inline member functions

   /*
   * Return q-field at beginning of block.
   */
   template <int D, class T> inline
   typename T::RField const& Propagator<D,T>::head() const
   {
      UTIL_CHECK(PropagatorTmplT::isSolved());
      return qFields_[0];
   }

   /*
   * Return q-field at end of block.
   */
   template <int D, class T> inline
   typename T::RField const& Propagator<D,T>::tail() const
   {
      UTIL_CHECK(PropagatorTmplT::isSolved());
      UTIL_CHECK(PolymerModel::isThread() || !PropagatorTmplT::isTailEnd());
      return qFields_[ns_-1];
   }

   /*
   * Return q-field at specified step.
   */
   template <int D, class T> inline
   typename T::RField const& Propagator<D,T>::q(int i) const
   {
      UTIL_CHECK(PropagatorTmplT::isSolved());
      return qFields_[i];
   }

   /*
   * Get the associated Block object by const reference.
   */
   template <int D, class T> inline
   typename T::Block const & Propagator<D,T>::block() const
   {
      UTIL_ASSERT(blockPtr_);
      return *blockPtr_;
   }

   /*
   * Get the associated Block object by non-const reference.
   */
   template <int D, class T> inline
   typename T::Block& Propagator<D,T>::block()
   {
      UTIL_ASSERT(blockPtr_);
      return *blockPtr_;
   }

   /*
   * Get the associated Mesh object by const reference.
   */
   template <int D, class T> inline
   Mesh<D> const & Propagator<D,T>::mesh() const
   {
      UTIL_ASSERT(meshPtr_);
      return *meshPtr_;
   }

   /*
   * Get the number of counter grid points.
   */
   template <int D, class T> inline
   int Propagator<D,T>::ns() const
   {  return ns_; }

   /*
   * Has memory been allocated for this propagator?
   */
   template <int D, class T> inline
   bool Propagator<D,T>::isAllocated() const
   {  return isAllocated_; }

   /*
   * Associate this propagator with a unique block.
   */
   template <int D, class T> inline
   void Propagator<D,T>::setBlock(typename T::Block& block)
   {  blockPtr_ = &block; }

} // namespace Rp
} // namespace Pscf
#endif
