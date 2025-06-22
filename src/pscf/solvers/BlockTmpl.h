#ifndef PSCF_BLOCK_TMPL_H
#define PSCF_BLOCK_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/chem/Edge.h>         // base class
#include <util/containers/Pair.h>   // member template

#include <cmath>

namespace Pscf
{ 

   using namespace Util;

   /**
   * Class template for a block solver in a block copolymer.
   *
   * Class template argument QT is a concrete propagator class. A 
   * BlockTmpl<QT> object has:
   *
   *   - two QT propagator objects, one per direction
   *   - a monomer concentration field for the block
   *   - a kuhn length
   * 
   * Each implementation of SCFT and/or FTS is defined in a different
   * enclosed namespace of Pscf. Each such implementation defines a
   * concrete propagator class and a concrete block class. By convention,
   * these are named Propagator and Block, respectively. The Block class 
   * in each implementation is derived from BlockTmpl<Propagator>, using 
   * the following syntax:
   * \code
   *  
   *    class Block : public BlockTmpl<Propagator>
   *    {
   *      ....
   *    }
   *
   * \endcode
   *
   * Design notes:
   *
   * The Block class of a polymer field theory implementation has member
   * variables that specify all the data that is needed to describe the 
   * block. This includes but is not limited to the monomer type id and
   * block length (defined by the Edge base class) and the kuhn length
   * (defined by this class template). In addition, each such Block class
   * may define a variety of member variables that define the numerical 
   * discretization of the block and any related parameters needed by 
   * the algorithm that solves the modified diffusion equation. 
   *   
   * The Block class will normally define one or more public functions
   * that can be called repeatedly by the Propagator::solve() function in
   * order to implement individual steps of the stepping algorithm used to
   * solve the MDE. 
   *  
   * The Block class also normally provides a void function named 
   * computeConcentration() that integrates the product of the two
   * associated propagators with respect to contour length in order
   * to compute the monomer concentration field associated with the
   * block.  This is normally called within the implementation of
   * solve() function of the associated Polymer class, within a loop 
   * over all blocks of the molecule that is called after solution of 
   * the modified diffusion equation for all propagators.
   *  
   * Here is an example of a simple interface of the Block class for an
   * implementation that is designed for self-consistent field theory:
   * \code
   *
   *   
   *   // -------------------------------------------------------------
   *   // One step of integration of the modified diffusion equation.
   *   //
   *   // \param in   input q-field from previous step
   *   // \param out  output q-field at next step
   *   // -------------------------------------------------------------
   *   void step(Propagator::QField const & in, Propagator::QField& out);
   *
   *   // -------------------------------------------------------------
   *   // Compute monomer concentration field for this block.
   *   // 
   *   // \param prefactor  numerical prefactor of phi/(Q*length)
   *   // -------------------------------------------------------------
   *   void computeConcentration(double prefactor);
   *
   * \endcode
   *
   * The step and computeConcentration functions in this example can both 
   * use private variables that depend on the monomer type and contour 
   * length of a particular block, and that apply to both of the two 
   * associated propagators.  Parameters that depend upon the step size 
   * used to discretize the MDE for a particular block generally cannot, 
   * however, be re-used by other blocks, because the MDE solver algorithm 
   * for a thread model may use slightly different step sizes in different 
   * blocks in order to divide each block into an even integer number of 
   * steps of equal length. 
   * 
   * \ingroup Pscf_Solver_Module
   */
   template <class QT>
   class BlockTmpl : public Edge
   {

   public:

      /**
      * Modified diffusion equation solver (propagator) type.
      */
      typedef QT PropagatorT;

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
      QT& propagator(int directionId);
   
      /**
      * Get a const Propagator for a specified direction.
      *
      * See above for number conventions.
      *
      * \param directionId integer index for direction (0 or 1)
      */
      QT const & propagator(int directionId) const;
   
      /**
      * Get the associated monomer concentration field.
      */
      typename QT::CFieldT& cField();

      /**
      * Get the associated const monomer concentration field.
      */
      typename QT::CFieldT const & cField() const;
   
      /**
      * Get monomer statistical segment length.
      */
      double kuhn() const;
  
   private:

      /// Pair of Propagator objects (one for each direction).
      Pair<PropagatorT> propagators_;

      /// Monomer concentration field.
      typename QT::CFieldT cField_;

      /// Monomer statistical segment length.
      double kuhn_;

   };

   // Inline member functions

   /*
   * Get a Propagator indexed by direction.
   */
   template <class QT>
   inline 
   QT& BlockTmpl<QT>::propagator(int directionId)
   {  return propagators_[directionId]; }

   /*
   * Get a const Propagator indexed by direction.
   */
   template <class QT>
   inline 
   QT const & BlockTmpl<QT>::propagator(int directionId) const
   {  return propagators_[directionId]; }

   /*
   * Get the monomer concentration field.
   */
   template <class QT>
   inline
   typename QT::CFieldT& BlockTmpl<QT>::cField()
   {  return cField_; }

   /*
   * Get the const monomer concentration field.
   */
   template <class QT>
   inline
   typename QT::CFieldT const & BlockTmpl<QT>::cField() const
   {  return cField_; }

   /*
   * Get the monomer statistical segment length. 
   */
   template <class QT>
   inline double BlockTmpl<QT>::kuhn() const
   {  return kuhn_; }

   // Non-inline functions

   /*
   * Constructor.
   */
   template <class QT>
   BlockTmpl<QT>::BlockTmpl()
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
   template <class QT>
   BlockTmpl<QT>::~BlockTmpl()
   {}

   /*
   * Set the monomer statistical segment length.
   */
   template <class QT>
   void BlockTmpl<QT>::setKuhn(double kuhn)
   {  kuhn_ = kuhn; }

}
#endif
