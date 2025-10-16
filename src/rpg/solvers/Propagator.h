#ifndef RPG_PROPAGATOR_H
#define RPG_PROPAGATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/solvers/PropagatorTmpl.h> // base class template

#include <prdc/cuda/RField.h>            // member
#include <pscf/cuda/DeviceArray.h>       // member 
#include <util/containers/DArray.h>      // member array


// Forward declarations
namespace Pscf { 
   template <int D> class Mesh; 
   namespace Rpg {
      template <int D> class Block;
   }
}

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * MDE solver for one-direction of one block.
   *
   * A fully initialized Propagator<D> has an association with a Block<D>
   * object that owns this propagator and its partner, and with a partner
   * Propagator<D> that solves the MDE within the same block in the
   * opposite direction. It also has an association with a Mesh<D> that 
   * describes a spatial grid, and associations with zero or more source
   * Propagator<D> objects that are used to compute an initial condition
   * for this propagator at the head vertex.
   *
   * The associated Block<D> stores information required to numerically
   * solve the modified diffusion equation (MDE), including quantities
   * that depend on the w-field associated with this block, the unit
   * cell parameters and (in the thread model) the contour step size.
   * These quantities are set and stored by the block because their values 
   * are the same for the two propagators owned by each block, but may be
   * different for different blocks. The algorithm used by a Propagator
   * to solve the the MDE repeatedly calls the step functions provided 
   * by the parent Block. 
   *
   * \ingroup Rpg_Solver_Module
   */
   template <int D>
   class Propagator : public PropagatorTmpl< Propagator<D> >
   {

   public:

      // Public typename aliases 

      /**
      * Base class type (partial template specialization).
      */
      using Base = PropagatorTmpl< Propagator<D> >;
      

      /**
      * Field type (function of position, defined on real space grid).
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
      * Allocate propagator arrays.
      * 
      * The parameter ns is the number of values of s at which q(r,s) is
      * calculated, including end values at the terminating vertices. 
      * See docs for the function ns(), which returns this value.
      *
      * The address of the associated Mesh<D> object is retained.
      *
      * An Exception is thrown if the propagator is already allocated.
      *
      * \param ns number of points along chain contour
      * \param mesh spatial discretization mesh
      */
      void allocate(int ns, const Mesh<D>& mesh);

      /**
      * Reallocate memory used by this propagator.
      *
      * This function is used when the value of ns is changed, which can
      * occur during some parameter sweeps. See docs for allocate and ns.
      *
      * An Exceptions is thrown if the propagator has not been previously
      * allocated.
      *
      * \param ns number of slices (including end points at vertices)
      */
      void reallocate(int ns);

      /**
      * Solve the modified diffusion equation (MDE) for this block.
      *
      * This function computes an initial q-field at the head of this
      * propagator, and then solves the modified diffusion equation for
      * the block to propagate from the head to the tail. Algorithms for
      * the thread or bead model may be used, depending the value of 
      * PolymerModel::model().
      */
      void solve();

      /**
      * Solve the MDE for a specified initial condition.
      *
      * This function solves the modified diffusion equation for this
      * block with a specified initial condition, which is given by the
      * function parameter "head". 
      *
      * \param head initial condition of q-field at head of block
      */
      void solve(RField<D> const & head);

      /**
      * Compute and return partition function for the molecule.
      *
      * This function computes the partition function Q for the molecule
      * as a spatial average of the pointwise product of the initial/head 
      * slice for this propagator and final/tail slice of its partner.
      */
      double computeQ();

      /**
      * Return q-field at specified slice.
      *
      * \param i step index
      */
      RField<D> const & q(int i) const;

      /**
      * Return q-field at initial (head) vertex.
      */
      RField<D> const & head();

      /**
      * Return q-field at terminal (tail) vertex.
      *
      * This function throws an Exception if invoked while the bead model 
      * is in use (i.e., if PolymerModel::isThread(() == false) and the tail
      * for this propagator is a chain end (i.e., if isTailEnd() == true)
      * In this case, tail slice is not needed, and so is not computed.
      */
      RField<D> const & tail() const;

      /**
      * Return the full array of q-fields as an unrolled 1D array.
      */
      DeviceArray<cudaReal> const & qAll();

      /**
      * Get the associated Block object by reference.
      */
      Block<D>& block();

      /**
      * Get the associated Block object by const reference.
      */
      Block<D> const & block() const;

      /**
      * Get the number of values of s (or slices), including head and tail.
      *
      * The value of ns is the number of values of s at which q(r,s) is
      * calculated, including the ned values at the terminating vertices
      * (the head and tail). In the bead model, this is two more than the
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

      /**
      * Array containing the entire propagator, stored on the device.
      *
      * The propagator data is stored contiguously to allow batched FFTs
      * to be performed on all contour steps simultaneously, which occurs
      * in Block::computeStress.
      */
      DeviceArray<cudaReal> qFieldsAll_;

      /**
      * Array of RFields, each associated with a slice of qFieldsAll_.
      *
      * Each RField<D> acts as a reference array that point to a slice
      * of qFieldsAll_. Each slice will have Mesh.size() elements and
      * corresponds to a single contour value. There are ns_ of such 
      * slices, which is the capacity of the outer DArray.
      */
      DArray< RField<D> > qFields_;

      /// Pointer to the associated Block.
      Block<D>* blockPtr_;

      /// Pointer to the associated Mesh
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

   };

   // Inline member functions

   /*
   * Return q-field at initial (head) vertex of block.
   */
   template <int D>
   inline
   RField<D> const & Propagator<D>::head()
   {
      UTIL_CHECK(isSolved());
      return qFields_[0];
   }

   /*
   * Return q-field at final (tail) vertex of block.
   */
   template <int D>
   inline
   RField<D> const & Propagator<D>::tail() const
   {
      UTIL_CHECK(isSolved());
      UTIL_CHECK(PolymerModel::isThread() || !isTailEnd());
      return qFields_[ns_-1];
   }

   /*
   * Return const q-field at specified step by reference.
   */
   template <int D>
   inline
   RField<D> const & Propagator<D>::q(int i) const
   {
      UTIL_CHECK(isSolved());
      return qFields_[i];
   }

   /*
   * Return the full array of q-fields.
   */
   template <int D>
   inline
   DeviceArray<cudaReal> const & Propagator<D>::qAll()
   {
      UTIL_CHECK(isSolved());
      return qFieldsAll_;
   }

   /*
   * Get the associated Block object (non-const reference)
   */
   template <int D>
   inline
   Block<D>& Propagator<D>::block()
   {
      assert(blockPtr_);
      return *blockPtr_;
   }

   /*
   * Get the associated Block object (const reference)
   */
   template <int D>
   inline
   Block<D> const & Propagator<D>::block() const
   {
      assert(blockPtr_);
      return *blockPtr_;
   }

   /*
   * Get the number ns of chain contour points.
   */
   template <int D>
   inline
   int Propagator<D>::ns() const
   {  return ns_; }

   template <int D>
   inline
   bool Propagator<D>::isAllocated() const
   {  return isAllocated_; }

   /*
   * Associate this propagator with a block and direction
   */
   template <int D>
   inline
   void Propagator<D>::setBlock(Block<D>& block)
   {  blockPtr_ = &block; }

   // Explicit instantiation declarations
   extern template class Propagator<1>;
   extern template class Propagator<2>;
   extern template class Propagator<3>;

}
}
#endif
