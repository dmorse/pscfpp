#ifndef PSCF_BLOCK_TMPL_H
#define PSCF_BLOCK_TMPL_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
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
   * Class TP is a concrete propagator class.
   *
   * \ingroup Pscf_Solvers_Module
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
      * Get the associated monomer concentration field.
      *
      * \param blockId integer index of associated block
      */
      typename TP::CField& cField();
   
      /**
      * Get monomer statistical segment length.
      */
      double kuhn() const;
  
   private:

      /// Propagators (one for each direction).
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
   * Get the monomer concentration field.
   */
   template <class TP>
   inline
   typename TP::CField& BlockTmpl<TP>::cField()
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
