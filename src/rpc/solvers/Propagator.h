#ifndef RPC_PROPAGATOR_H
#define RPC_PROPAGATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/solvers/PropagatorTmpl.h> // base class template
#include <prdc/cpu/RField.h>             // member template and typedef
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
   * A fully initialized Propagator<D> has associations with a 
   * Block<D> that owns this propagator and with a partner Propagator<D>.
   * It also has an association with a Mesh<D> that describes a spatial 
   * grid, and to source Propagator<D> objects that are used to compute
   * an initial condition for the propagator at the head vertex.
   *
   * The associated Block<D> stores information required to numerically
   * solve the modified diffusion equation (MDE) within the block, 
   * including the contour step size ds (in the thread model) and all 
   * parameters that depend on ds, unit cell parameters and the w-field
   * associated with this block. These quantities are set and stored by 
   * the block because their values must be the same for both of the two 
   * propagators owned by each block (i.e., this propagator and its 
   * partner). The algorithm used by a propagator to solve the MDE 
   * repeatedly calls step functions provided by the associated block.
   *
   * \ingroup Rpc_Solver_Module
   */
   template <int D>
   class Propagator : public PropagatorTmpl< Propagator<D> >
   {

   public:

      // Public typedefs

      /**
      * Generic field (function of position, defined on regular grid).
      */
      typedef RField<D> FieldT;

      /**
      * Chemical potential field type (r-grid format)
      */
      typedef RField<D> WFieldT;

      /**
      * Monomer concentration field type (r-grid format)
      */
      typedef RField<D> CFieldT;

      /**
      * Propagator q-field type, i.e., q(r,s) at fixed s.
      */
      typedef RField<D> QFieldT;

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
      * \param ns  number of slices (including end points)
      * \param mesh  spatial discretization mesh
      */ 
      void allocate(int ns, const Mesh<D>& mesh);

      /**
      * Reallocate memory used by this propagator.
      * 
      * This function is used when the value of ns is changed after initial
      * allocation. This occurs during parameter sweeps that change the
      * block length. See the docs for the function ns() for the definition
      * of ns.
      *
      * The spatial mesh is set by derefencing a pointer to the associated
      * Mesh<D> object, which was set by a previous call to allocate.
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
      * the function parameter "head". Algorithms for the thread or bead
      * model may be used, depending on value of PolymerModel::model().
      *
      * \param head  initial condition of q-field at head of block
      */
      void solve(QFieldT const & head);
 
      /**
      * Compute and return partition function for the polymer.
      *
      * This function computes the partition function Q for the 
      * molecule as a spatial average of pointwise product of the
      * initial/head Qfield for this propagator and the final/tail 
      * Qfield of its partner. An algorithm appropriate to the 
      * thread or bead model is chosen based on the value of
      * PolymerModel::model().
      *
      * \return value of Q (spatial average of q*q^{+} at head)
      */ 
      double computeQ();

      /**
      * Return q-field at specified step.
      *
      * \param i step index, 0 <= i < ns
      */
      const QFieldT& q(int i) const;

      /**
      * Return q-field at beginning of the block (initial condition).
      */
      const QFieldT& head() const;

      /**
      * Return q-field at the end of the block.
      */
      const QFieldT& tail() const;

      /**
      * Get the associated Block object by reference.
      */
      Block<D>& block();

      /**
      * Number of values of s (or slices), including head and tail.
      *
      * The value of ns is the number of values of s at which q(r,s) is
      * calculated, including the end values at the terminating vertices
      * (the head and tail).  This is one more than the number of 
      * contour variable steps. If the propagator does not own both
      * terminal vertex beads, it ns may exceed nBead.
      */
      int ns() const;

      /**
      * Has memory been allocated for this propagator?
      */
      bool isAllocated() const;

      // Inherited public members with non-dependent names

      using PropagatorTmpl< Propagator<D> >::nSource;
      using PropagatorTmpl< Propagator<D> >::source;
      using PropagatorTmpl< Propagator<D> >::partner;
      using PropagatorTmpl< Propagator<D> >::setIsSolved;
      using PropagatorTmpl< Propagator<D> >::isSolved;
      using PropagatorTmpl< Propagator<D> >::hasPartner;
      using PropagatorTmpl< Propagator<D> >::isHeadEnd;
      using PropagatorTmpl< Propagator<D> >::isTailEnd;

   private:
     
      /// Array of propagator slices at different contour variable values.
      DArray<QFieldT> qFields_;

      /// Workspace
      QFieldT work_;

      /// Pointer to associated Block.
      Block<D>* blockPtr_;

      /// Pointer to associated Mesh
      Mesh<D> const * meshPtr_;

      /// Number of grid points = # of contour length steps + 1
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
      void assign(QFieldT& lhs, QFieldT const & rhs);

   };

   // Inline member functions

   /*
   * Return q-field at beginning of block.
   */
   template <int D>
   inline 
   typename Propagator<D>::QFieldT const& Propagator<D>::head() const
   {  
      UTIL_CHECK(isSolved()); 
      return qFields_[0]; }

   /*
   * Return q-field at end of block.
   */
   template <int D>
   inline 
   typename Propagator<D>::QFieldT const& Propagator<D>::tail() const
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
   typename Propagator<D>::QFieldT const& Propagator<D>::q(int i) const
   {  
      UTIL_CHECK(isSolved()); 
      return qFields_[i]; 
   }

   /*
   * Get the associated Block object.
   */
   template <int D>
   inline 
   Block<D>& Propagator<D>::block()
   {
      assert(blockPtr_);  
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
      assert(blockPtr_);  
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
