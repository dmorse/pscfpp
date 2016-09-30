#ifndef FD1D_PROPAGATOR_H
#define FD1D_PROPAGATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/solvers/PropagatorTmpl.h> // base class template
#include <util/containers/DArray.h>      // member template

namespace Pscf { 
namespace Fd1d
{ 

   class Block;
   using namespace Util;

   /**
   * MDE solver for one-direction of one block.
   *
   * \ingroup Pscf_Fd1d_Module
   */
   class Propagator : public PropagatorTmpl<Propagator>
   {

   public:

      // Public typedefs

      /**
      * Generic field (function of position).
      */
      typedef DArray<double> Field;

      /**
      * Chemical potential field type.
      */
      typedef DArray<double> WField;

      /**
      * Monomer concentration field type.
      */
      typedef DArray<double> CField;

      /**
      * Propagator q-field type.
      */
      typedef DArray<double> QField;

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
      void setBlock(Block& block);

      /**
      * Associate this propagator with a block.
      * 
      * \param ns number of contour length steps
      * \param nx number of spatial steps
      */ 
      void allocate(int ns, int nx);

      /**
      * Solve the modified diffusion equation (MDE) for this block.
      *
      * This function computes an initial QField at the head of this
      * block, and then solves the modified diffusion equation for 
      * the block to propagate from the head to the tail. The initial
      * QField at the head is computed by pointwise multiplication of
      * of the tail QFields of all source propagators.
      */
      void solve();
  
      /**
      * Solve the MDE for a specified initial condition.
      *
      * This function solves the modified diffusion equation for this
      * block with a specified initial condition, which is given by 
      * head parameter of the function. The function is intended for 
      * use in testing.
      *
      * \param head initial condition of QField at head of block
      */
      void solve(const QField& head);
 
      /**
      * Compute and return partition function for the molecule.
      *
      * This function computes the partition function Q for the 
      * molecule as a spatial average of the initial/head Qfield 
      * for this propagator and the final/tail Qfield of its
      * partner. 
      */ 
      double computeQ();

      /**
      * Return q-field at specified step.
      *
      * \param i step index
      */
      const QField& q(int i) const;

      /**
      * Return q-field at beginning of block (initial condition).
      */
      const QField& head() const;

      /**
      * Return q-field at end of block.
      */
      const QField& tail() const;

      /**
      * Get the associated Block object by reference.
      */
      Block & block();

      /**
      * Has memory been allocated for this propagator?
      */
      bool isAllocated() const;

   protected:

      /**
      * Compute initial QField at head from tail QFields of sources.
      */
      void computeHead();

   private:
     
      // Array of statistical weight fields 
      DArray<QField> qFields_;

      // Workspace
      QField work_;

      /// Pointer to associated Block.
      Block* blockPtr_;

      /// Number of contour length steps = # grid points - 1.
      int ns_;

      /// Number of spatial grid points.
      int nx_;

      /// Is this propagator allocated?
      bool isAllocated_;

   };

   // Inline member functions

   /*
   * Return q-field at beginning of block.
   */
   inline Propagator::QField const& Propagator::head() const
   {  return qFields_[0]; }

   /*
   * Return q-field at end of block, after solution.
   */
   inline Propagator::QField const& Propagator::tail() const
   {  return qFields_[ns_-1]; }

   /*
   * Return q-field at specified step.
   */
   inline Propagator::QField const& Propagator::q(int i) const
   {  return qFields_[i]; }

   /*
   * Get the associated Block object.
   */
   inline Block& Propagator::block()
   {
      assert(blockPtr_);  
      return *blockPtr_; 
   }

   /*
   * Associate this propagator with a block and direction
   */
   inline void Propagator::setBlock(Block& block)
   {  blockPtr_ = &block; }

}
} 
#endif
