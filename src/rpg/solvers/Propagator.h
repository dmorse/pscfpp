#ifndef RPG_PROPAGATOR_H
#define RPG_PROPAGATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <prdc/cuda/RField.h>            // member array
#include <pscf/cuda/DeviceArray.h>       // member array
#include <pscf/solvers/PropagatorTmpl.h> // base class template
#include <util/containers/DArray.h>      // member array

namespace Pscf { 

   // Forward declaration
   template <int D> class Mesh;

namespace Rpg { 

   // Forward declaration
   template <int D> class Block;

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cuda;

   /**
   * MDE solver for one-direction of one block.
   *
   * A fully initialized Propagator<D> has an association with a 
   * Block<D> that owns this propagator and its partner, and has an 
   * association with a Mesh<D> that describes a spatial grid, in 
   * addition to associations with partner and source Propagator<D>
   * objects that are managed by the PropagatorTmpl base class template. 
   *
   * The associated Block<D> stores information required to numerically
   * solve the modified diffusion equation (MDE), including the contour
   * step size ds and all parameters that depend on ds. These quantities 
   * are set and stored by the block because their values must be the 
   * same for the two propagators owned by each block (i.e., this 
   * propagator and its partner). The algorithm used by a propagator 
   * to solve the the MDE simply repeatedly calls the step() function 
   * of the associated block, because that function has access to all 
   * the parameters used in the numerical solution.
   *
   * \ingroup Rpg_Solvers_Module
   */
   template <int D>
   class Propagator : public PropagatorTmpl< Propagator<D> >
   {

   public:

      // Public typedefs (used by template classes)

      /**
      * Generic field (function of position).
      */ 
      typedef RField<D> Field;

      /**
      * Chemical potential field type.
      */ 
      typedef RField<D> WField;

      /**
      * Monomer concentration field type.
      */
      typedef RField<D> CField;

      /**
      * Propagator q-field type.
      */
      typedef RField<D> QField;

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
      * \param ns number of points along chain contour
      * \param mesh spatial discretization mesh
      */ 
      void allocate(int ns, const Mesh<D>& mesh);

      /**
      * Reallocate memory used by this propagator.
      * 
      * This function is used when the value of ns is changed,
      * which occurs during some parameter sweeps. 
      * 
      * The parameter ns is the number of values of s at which
      * q(r,s) is calculated, including the end values at the
      * terminating vertices. This is one more than the number 
      * of contour variable steps. 
      * 
      * \param ns number of points along chain contour
      */ 
      void reallocate(int ns);

      /**
      * Solve the modified diffusion equation (MDE) for this block.
      *
      * This function computes an initial QField at the head of this
      * propagator, and then solves the modified diffusion equation for 
      * the block to propagate from the head to the tail. The initial
      * QField at the head is computed by pointwise multiplication of
      * of the tail QFields of all source propagators.
      */
      void solve();
  
      /**
      * Solve the MDE for a specified initial condition.
      *
      * This function solves the modified diffusion equation for this
      * block with a specified initial condition, which is given by head
      * parameter of the function. The function is intended for use in
      * testing.
      *
      * \param head initial condition of QField at head of block
      */
      void solve(RField<D> const & head);
 
      /**
      * Compute and return partition function for the molecule.
      *
      * This function computes the partition function Q for the molecule
      * as a spatial average of the product of the initial / head Qfield 
      * for this propagator and the final / tail Qfield of its partner. 
      */ 
      double computeQ();
      
      /**
      * Return const q-field at specified step by reference (after solving).
      *
      * \param i step index
      */
      RField<D> const & q(int i) const;

      /**
      * Return q-field at beginning of block (initial condition).
      */
      RField<D> const & head();

      /**
      * Return q-field at end of block (after propagator is solved).
      */
      RField<D> const & tail() const;

      /**
      * Return the full array of q-fields (after propagator is solved).
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
      * Get the number of chain contour points.
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
      using PropagatorTmpl< Propagator<D> >::ownsHead;
      using PropagatorTmpl< Propagator<D> >::ownsTail;

   protected:

      /**
      * Compute initial QField at head for the thread model. 
      *
      * In the thread model, the head slice of each propagator is the
      * product of tail slices for incoming propagators from other bonds
      * that terminate at the head vertex.
      */
      void computeHeadThread();

      /**
      * Compute initial QField at head for the bead model.
      *
      * In the bond model, the head slice for each propagator is given
      * by the product of tails slices for source propagators, times an
      * additional bond weight exp(-W(r)*ds) if this propagator owns the
      * head vertex bead.
      */
      void computeHeadBead();

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
      * These RFields will act as reference arrays that point to slices of
      * qFieldsAll_. Each slice will have Mesh.size() elements and 
      * correspond to a single contour point, and there will be ns_ of
      * these reference arrays.
      */
      DArray<RField<D> > qFields_;

      /// Pointer to associated Block.
      Block<D>* blockPtr_;

      /// Pointer to associated Mesh
      Mesh<D> const * meshPtr_;

      /// Number of chain contour positions (= # contour steps + 1).
      int ns_;

      /// Is this propagator allocated?
      bool isAllocated_;

   };

   // Inline member functions

   /*
   * Return q-field at beginning of block.
   */
   template <int D>
   inline 
   RField<D> const & Propagator<D>::head()
   {  
      UTIL_CHECK(isAllocated());
      return qFields_[0];
   }

   /*
   * Return q-field at end of block, after solution.
   */
   template <int D>
   inline 
   RField<D> const & Propagator<D>::tail() const
   {  
      UTIL_CHECK(isSolved());
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
   * Get the associated Block object (non-const reference)
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


   #ifndef RPG_PROPAGATOR_TPP
   // Suppress implicit instantiation
   extern template class Propagator<1>;
   extern template class Propagator<2>;
   extern template class Propagator<3>;
   #endif

}
}
#endif
