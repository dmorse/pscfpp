#ifndef PSCF_BLOCK_TMPL_H
#define PSCF_BLOCK_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/BlockDescriptor.h> // base class
#include <util/containers/Pair.h>      // member template

#include <cmath>

namespace Pscf
{ 

   using namespace Util;

   /**
   * Class template for a block in a block copolymer.
   *
   * Class TP is a concrete propagator class. A BlockTmpl<TP> object
   * has:
   *
   *   - two TP propagator objects, one per direction
   *   - a single monomer concentration field 
   *   - a single kuhn length
   * 
   * Each implementation of self-consistent field theory (SCFT) is defined
   * in a different sub-namespace of Pscf. Each such implementation defines
   * a concrete propagator class and a concrete block class, which ar 
   * named Propagator and Block by convention.The Block class in each 
   * implementation is derived from BlockTmpl<Propagator>, using the 
   * following syntax:
   * \code
   *  
   *    class Block : public BlockTmpl<Propagator>
   *    {
   *      ....
   *    }
   *
   * \endcode
   * The algorithms for taking one step of integration of the modified
   * diffusion equation and for computing the monomer concentration field
   * arising from monomers in one block must be implemented in the Block
   * class. These algorithms must be implemented in member functions of
   * the concrete Block class with the following interfaces:
   * \code
   *
   *   
   *   // ---------------------------------------------------------------
   *   // Take one step of integration of the modified diffusion equation.
   *   // ---------------------------------------------------------------
   *   void step(Propagator::QField const & in, Propagator::QField& out);
   *
   *   // ---------------------------------------------------------------
   *   // Compute monomer concentration field for this block.
   *   // 
   *   // \param prefactor  numerical prefactor of phi/(Q*length)
   *   // ---------------------------------------------------------------
   *   void computeConcentration(double prefactor);
   *
   * \endcode
   * These core algorithms are implemented in the Block class, rather than 
   * the Propagator class, because the data required to implement these 
   * algorithms generally depends on the monomer type and contour length 
   * step size ds of a particular block, and can thus be shared by the two 
   * propagators associated with a particular block. The data required 
   * to implement these algorithms cannot, however, be shared among 
   * propagators in different blocks of the same monomer type because the 
   * requirement that the length of each block be divided into an integer 
   * number of contour length steps implies that different blocks of 
   * arbitrary user-specified length must generally be assumed to have 
   * slightly different values for the step size ds.
   * 
   * The step() function is called in the implementation of the 
   * PropagatorTmpl::solve() member function, within a loop over steps.
   * The computeConcentration() function is called in the implementation
   * of the PolymerTmpl::solve() member function, within a loop over all
   * blocks of the molecule that is called after solution of the modified
   * diffusion equation for all propagators.
   *
   * \ingroup Pscf_Solver_Module
   */
   template <class TP>
   class BlockTmpl : public BlockDescriptor
   {

   public:

      // Modified diffusion equation propagator for one block.
      typedef TP Propagator;

      // Monomer concentration field.
      typedef typename TP::CField CField;
 
      // Chemical potential field.
      typedef typename TP::WField WField;

      /**
      * Constructor.
      */
      BlockTmpl();
 
      /**
      * Destructor.
      */
      virtual ~BlockTmpl();

      /**
      * Set monomer statistical segment length.
      *
      * \param kuhn monomer statistical segment length
      */
      virtual void setKuhn(double kuhn);

      /**
      * Get a Propagator for a specified direction.
      *
      * For a block with v0 = vertexId(0) and v1 = vertexId(1),
      * propagator(0) propagates from vertex v0 to v1, while
      * propagator(1) propagates from vertex v1 to v0.
      *
      * \param directionId integer index for direction (0 or 1)
      */
      TP& propagator(int directionId);
   
      /**
      * Get a const Propagator for a specified direction.
      *
      * See above for number conventions.
      *
      * \param directionId integer index for direction (0 or 1)
      */
      TP const & propagator(int directionId) const;
   
      /**
      * Get the associated monomer concentration field.
      */
      typename TP::CField& cField();

      /**
      * Get the associated const monomer concentration field.
      */
      typename TP::CField const & cField() const;
   
      /**
      * Get monomer statistical segment length.
      */
      double kuhn() const;
  
   private:

      /// Pair of Propagator objects (one for each direction).
      Pair<Propagator> propagators_;

      /// Monomer concentration field.
      CField cField_;

      /// Monomer statistical segment length.
      double kuhn_;

   };

   // Inline member functions

   /*
   * Get a Propagator indexed by direction.
   */
   template <class TP>
   inline 
   TP& BlockTmpl<TP>::propagator(int directionId)
   {  return propagators_[directionId]; }

   /*
   * Get a const Propagator indexed by direction.
   */
   template <class TP>
   inline 
   TP const & BlockTmpl<TP>::propagator(int directionId) const
   {  return propagators_[directionId]; }

   /*
   * Get the monomer concentration field.
   */
   template <class TP>
   inline
   typename TP::CField& BlockTmpl<TP>::cField()
   {  return cField_; }

   /*
   * Get the const monomer concentration field.
   */
   template <class TP>
   inline
   typename TP::CField const & BlockTmpl<TP>::cField() const
   {  return cField_; }

   /*
   * Get the monomer statistical segment length. 
   */
   template <class TP>
   inline double BlockTmpl<TP>::kuhn() const
   {  return kuhn_; }

   // Non-inline functions

   /*
   * Constructor.
   */
   template <class TP>
   BlockTmpl<TP>::BlockTmpl()
    : propagators_(),
      cField_(),
      kuhn_(0.0)
   {
      propagator(0).setDirectionId(0);
      propagator(1).setDirectionId(1);
      propagator(0).setPartner(propagator(1));
      propagator(1).setPartner(propagator(0));
   }

   /*
   * Destructor.
   */
   template <class TP>
   BlockTmpl<TP>::~BlockTmpl()
   {}

   /*
   * Set the monomer statistical segment length.
   */
   template <class TP>
   void BlockTmpl<TP>::setKuhn(double kuhn)
   {  kuhn_ = kuhn; }

}
#endif
