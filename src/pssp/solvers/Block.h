#ifndef PSSP_BLOCK_H
#define PSSP_BLOCK_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.h"                   // base class argument
#include <pscf/solvers/BlockTmpl.h>       // base class template

namespace Pscf { 
namespace Pssp { 

   // class Domain;
   using namespace Util;

   /**
   * Block within a branched polymer.
   *
   * Derived from BlockTmpl<Propagator>. A BlockTmpl<Propagator> has two 
   * Propagator members and is derived from BlockDescriptor.
   *
   * \ingroup Pscf_Fd1d_Module
   */
   template <int D>
   class Block : public BlockTmpl< Propagator<D> >
   {

   public:

      /**
      * Generic field (base class)
      */
      typedef typename Propagator<D>::Field Field;

      /**
      * Monomer chemical potential field.
      */
      typedef typename Propagator<D>::WField WField;

      /**
      * Constrained partition function q(r,s) for fixed s.
      */
      typedef typename Propagator<D>::QField QField;

      // Member functions

      /**
      * Constructor.
      */
      Block();

      /**
      * Destructor.
      */
      ~Block();

      /**
      * Initialize discretization and allocate required memory.
      *
      * \param domain associated Domain object, with grid info
      * \param ds desired (optimal) value for contour length step
      */
      // void setDiscretization(Domain const & domain, double ds);
      void setDiscretization(double ds);

      /**
      * Set solver for this block.
      */
      void setupSolver(WField const & w);

      /**
      * Compute unnormalized concentration for block by integration.
      *
      * Upon return, grid point r of array cField() contains the 
      * integral int ds q(r,s)q^{*}(r,L-s) times the prefactor, 
      * where q(r,s) is the solution obtained from propagator(0), 
      * and q^{*} is the solution of propagator(1),  and s is
      * a contour variable that is integrated over the domain 
      * 0 < s < length(), where length() is the block length.
      *
      * \param prefactor multiplying integral
      */ 
      void computeConcentration(double prefactor);

      /**
      * Compute step of integration loop, from i to i+1.
      */
      void step(QField const & q, QField& qNew);

      #if 0
      /**
      * Return associated domain by reference.
      */
      Domain const & domain() const;
      #endif

      /**
      * Number of contour length steps.
      */
      int ns() const;

      using BlockTmpl< Propagator<D> >::setKuhn;
      using BlockTmpl< Propagator<D> >::setupSolver;
      using BlockTmpl< Propagator<D> >::computeConcentration;
      using BlockTmpl< Propagator<D> >::propagator;
      using BlockTmpl< Propagator<D> >::cField;
      using BlockTmpl< Propagator<D> >::length;
      using BlockTmpl< Propagator<D> >::kuhn;

      using BlockDescriptor::setId;
      using BlockDescriptor::setVertexIds;
      using BlockDescriptor::setMonomerId;
      using BlockDescriptor::setLength;
      using BlockDescriptor::id;
      using BlockDescriptor::monomerId;
      using BlockDescriptor::vertexIds;
      using BlockDescriptor::vertexId;
      using BlockDescriptor::length;

   private:

      /// Pointer to associated Domain object.
      // Domain const * domainPtr_;

      /// Contour length step size.
      double ds_;

      /// Number of contour length steps = # grid points - 1.
      int ns_;

   };

   // Inline member functions

   #if 0
   /// Get Domain by reference.
   template <int D>
   inline Domain const & Block<D>::domain() const
   {   
      UTIL_ASSERT(domainPtr_);
      return *domainPtr_;
   }
   #endif

   /// Get number of contour steps.
   template <int D>
   inline int Block<D>::ns() const
   {  return ns_; }

}
}
#include "Block.tpp"
#endif
