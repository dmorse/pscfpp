#ifndef RPC_PROPAGATOR_H
#define RPC_PROPAGATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/solvers/PropagatorTmpl.h> // base class template
#include <prdc/cpu/RField.h>             // member template 
#include <util/containers/DArray.h>      // member template
#include <util/containers/FArray.h>      // member template

namespace Pscf { template <int D> class Mesh; }

namespace Pscf { 
namespace Rpc { 

   // Forward declaration
   template <int D> class Block;

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /**
   * MDE solver for one direction of one block.
   *
   * A fully initialized Propagator<D> has an associations with a Block<D>
   * object that owns this propagator and its partner, and with a partner 
   * Propagator<D> that solves the MDE within the same block in the 
   * opposite direction. It also has an association with a Mesh<D> that 
   * describes a spatial grid, and associations with zero or more source 
   * Propagator<D> objects that are used to compute an initial condition 
   * for this propagator at the head vertex.
   *
   * The associated Block<D> stores information required to numerically
   * solve the modified diffusion equation (MDE), including quantities 
   * that depend upon the w-field associated with this block, the unit
   * cell parameters and (in the thread model) the contour step size.
   * These quantities are set and stored by the block because their values 
   * are the same for the two propagators owned by each block, but may be
   * different for different blocks.  The algorithm used by a Propagator 
   * to solve the MDE repeatedly calls step functions provided by the 
   * parent Block.
   *
   * \ingroup Rpc_Solver_Module
   */
   template <int D>
   class Propagator : public PropagatorTmpl< Propagator<D> >
   {

   public:

      // Public typename aliases

      /**
      * Base class (partial template specialization).
      */
      using Base = PropagatorTmpl< Propagator<D> >;

      /**
      * Field type (function of position, defined on a r-space grid).
      */
      using FieldT = RField<D>;

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
      * Associate this propagator with a block.
      *
      * \param block associated Block object.
      */ 
      void setBlock(Block<D>& block);

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
      void allocate(int ns, const Mesh<D>& mesh);

      /**
      * Reallocate memory used by this propagator.
      * 
      * This function is used when the value of ns is changed, which can
      * occur during some parameter sweeps. See docs for allocate and ns.
      * 
      * An Exception is thrown if the propagator has not been previously
      * allocated, or if the parameter ns is equal to the current value.
      * 
      * \param ns  number of slices (including end points)
      */ 
      void reallocate(int ns);

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
      void solve(FieldT const & head);
 
      /**
      * Compute and return partition function for the polymer.
      *
      * This function computes the partition function Q for the molecule
      * as a spatial average of the pointwise product of the initial/head 
      * slice for this propagator and the final/tail slice of its partner. 
      *
      * \return value of Q (spatial average of q*q^{+} at head)
      */ 
      double computeQ() const;

      /**
      * Return q-field at specified step.
      *
      * \param i step index, 0 <= i < ns
      */
      const FieldT& q(int i) const;

      /**
      * Return q-field at beginning of the block (initial condition).
      */
      const FieldT& head() const;

      /**
      * Return q-field at the end of the block.
      *
      * This function throws an Exception if invoked while the bead model
      * is in use (i.e., if PolymerModel::isThread() == false) and the tail 
      * for this propagator is a chain end (i.e., if isTailEnd() == true).
      * In this case, the tail slice is not needed, and so is not computed.
      */
      const FieldT& tail() const;

      /**
      * Get the associated Block object by const reference.
      */
      Block<D> const & block() const;

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

      // Inherited public members with non-dependent names

      using Base::nSource;
      using Base::source;
      using Base::partner;
      using Base::setIsSolved;
      using Base::isSolved;
      using Base::hasPartner;
      using Base::isHeadEnd;
      using Base::isTailEnd;

   private:
     
      /// Array of propagator slices at different contour variable values.
      DArray<FieldT> qFields_;

      /// Workspace.
      FieldT work_;

      /// Pointer to the associated Block.
      Block<D>* blockPtr_;

      /// Pointer to the associated Mesh.
      Mesh<D> const * meshPtr_;

      /// Number of slices, including head and tail slices.
      int ns_;

      /// Is this propagator allocated?
      bool isAllocated_;

      /**
      * Compute initial q-field at head.
      * 
      * In either model, the head slice of each propagator is the product
      * of tail slices for incoming propagators from other bonds that
      * terminate at the head vertex.
      */
      void computeHead();

      /**
      * Assign one slice to another (RHS = LHS).
      */
      void assign(FieldT& lhs, FieldT const & rhs);

   };

   // Inline member functions

   /*
   * Return q-field at beginning of block.
   */
   template <int D>
   inline 
   typename Propagator<D>::FieldT const& Propagator<D>::head() const
   {  
      UTIL_CHECK(isSolved()); 
      return qFields_[0]; }

   /*
   * Return q-field at end of block.
   */
   template <int D>
   inline 
   typename Propagator<D>::FieldT const& Propagator<D>::tail() const
   {
      UTIL_CHECK(isSolved()); 
      UTIL_CHECK(PolymerModel::isThread() || !isTailEnd());
      return qFields_[ns_-1]; 
   }

   /*
   * Return q-field at specified step.
   */
   template <int D>
   inline 
   typename Propagator<D>::FieldT const& Propagator<D>::q(int i) const
   {  
      UTIL_CHECK(isSolved()); 
      return qFields_[i]; 
   }

   /*
   * Get the associated Block object by const reference.
   */
   template <int D>
   inline 
   Block<D> const & Propagator<D>::block() const
   {
      UTIL_ASSERT(blockPtr_);  
      return *blockPtr_; 
   }

   /*
   * Get the number of counter grid points.
   */
   template <int D>
   inline int Propagator<D>::ns() const
   {  return ns_; }

   template <int D>
   inline bool Propagator<D>::isAllocated() const
   {  return isAllocated_; }

   /*
   * Associate this propagator with a unique block.
   */
   template <int D>
   inline void Propagator<D>::setBlock(Block<D>& block)
   {
      UTIL_ASSERT(blockPtr_);  
      blockPtr_ = &block; 
   }

   #ifndef RPC_PROPAGATOR_TPP
   // Suppress implicit instantiation
   extern template class Propagator<1>;
   extern template class Propagator<2>;
   extern template class Propagator<3>;
   #endif

}
}
#endif
