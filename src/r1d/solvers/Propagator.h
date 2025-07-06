#ifndef R1D_PROPAGATOR_H
#define R1D_PROPAGATOR_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/solvers/PropagatorTmpl.h> // base class template
#include <util/containers/DArray.h>      // member template

namespace Pscf { 
namespace R1d
{ 

   class Block;
   using namespace Util;

   /**
   * MDE solver for one-direction of one block.
   *
   * \ingroup R1d_Solver_Module
   */
   class Propagator : public PropagatorTmpl<Propagator>
   {

   public:

      // Public typedefs

      /**
      * Generic field (function of position).
      */
      typedef DArray<double> FieldT;

      /**
      * Propagator q-field type.
      */
      typedef DArray<double> QFieldT;

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
      * Set discretization and allocate memory.
      * 
      * \param ns number of contour length steps
      * \param nx number of spatial steps
      */ 
      void allocate(int ns, int nx);

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
      * \param ns number of slices (including end points)
      */ 
      void reallocate(int ns);

      /**
      * Solve the modified diffusion equation (MDE) for this block.
      *
      * This function computes an initial q-field at the head of this
      * block, and then solves the modified diffusion equation for 
      * the block to propagate from the head to the tail. The initial
      * q-field at the head is computed by pointwise multiplication of
      * of the tail q-field of all source propagators.
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
      * \param head initial condition of q-field at head of block
      */
      void solve(const QFieldT& head);
 
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
      QFieldT const & q(int i) const;

      /**
      * Return q-field at beginning of block (initial condition).
      */
      QFieldT const & head() const;

      /**
      * Return q-field at end of block.
      */
      QFieldT const & tail() const;
      
      /**
      * Number of values of s (or slices), including head and tail.
      *
      * The value of ns is the number of values of s at which q(r,s) is
      * calculated, including the end values at the terminating vertices
      * (the head and tail).  This is one more than the number of contour 
      * variable steps. 
      */
      int ns() const;

      /**
      * Has memory been allocated for this propagator?
      */
      bool isAllocated() const;

   private:
     
      /// Array of statistical weight fields 
      DArray<QFieldT> qFields_;

      /// Workspace
      QFieldT work_;

      /// Pointer to associated Block.
      Block* blockPtr_;

      /// Number of contour length steps = # grid points - 1.
      int ns_;

      /// Number of spatial grid points.
      int nx_;

      /// Is this propagator allocated?
      bool isAllocated_;

      /**
      * Get the associated Block object by reference.
      */
      Block & block();

      /**
      * Compute initial q-field at head from tail q-fields of sources.
      */
      void computeHead();

   };

   // Inline member functions

   /*
   * Return q-field at beginning of block.
   */
   inline Propagator::QFieldT const& Propagator::head() const
   {  return qFields_[0]; }

   /*
   * Return q-field at end of block, after solution.
   */
   inline Propagator::QFieldT const& Propagator::tail() const
   {  return qFields_[ns_-1]; }

   /*
   * Return q-field at specified step.
   */
   inline Propagator::QFieldT const& Propagator::q(int i) const
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
   * Get the number of counter grid points.
   */
   inline int Propagator::ns() const
   {  return ns_; }

   /*
   * Associate this propagator with a block and direction
   */
   inline void Propagator::setBlock(Block& block)
   {  blockPtr_ = &block; }

}
} 
#endif
