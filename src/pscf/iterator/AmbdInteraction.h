#ifndef PSCF_AMBD_INTERACTION_H
#define PSCF_AMBD_INTERACTION_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/inter/Interaction.h>       // base class
#include <util/containers/Array.h>        // argument (template)
#include <util/containers/Matrix.h>       // argument (template)
#include <util/global.h>                  

namespace Pscf {

   using namespace Util;

   /**
   * Modified interaction to compute residual defn. of Arora et al.
   *
   * This class computes and provides access to copies of chi and as
   * well as several other auxiliary quantities that are needed to 
   * compute the SCFT residual for multi-component systems as defined
   * by Arora, Morse, Bates and Dorfman (AMBD) (JCP 2017) for use
   * in Anderson-mixing iteration (see reference below). The class
   * is designed to be used as a augmented local copy of the
   * interaction class within iterator classes that use this residual
   * definition. When used in this way, it should be updated just 
   * before entering the iteration loop.
   * 
   * Reference:
   * A. Arora, D.C. Morse, F.S. Bates and K.D Dorfman,
   * J. Chem. Phys vol. 146, 244902 (2017).
   * 
   * \ingroup Pscf_Iterator_Module
   */
   class AmbdInteraction 
   {

   public:

      /**
      * Constructor.
      */
      AmbdInteraction();

      /**
      * Destructor.
      */
      virtual ~AmbdInteraction();

      /**
      * Set number of monomers and allocate required memory.
      *
      * \param nMonomer number of different monomer types
      */
      void setNMonomer(int nMonomer);

      /**
      * Update all computed quantities.
      *
      * \param interaction Interaction object with current chi matrix
      */
      void update(Interaction const & interaction);

      /**
      * Return one element of the chi matrix.
      *
      * \param i row index
      * \param j column index
      */
      double chi(int i, int j) const;

      /**
      * Return one element of the inverse chi matrix.
      *
      * \param i row index
      * \param j column index
      */
      double chiInverse(int i, int j) const;

      /** 
      * Return one element of the potent matrix P.
      *
      * This matrix is defined by AMBD as:
      * 
      *     P = I - e e^{T} chi^{-1} /S
      *
      * where I is the Nmonomer x Nmonomer identity matrix, the column
      * vector e is a vector for which e_i = 1 for all i, and S is
      * the sum of the N^{2} elements of the matrix chi^{-1}, also 
      * given by S = e^{T} chi^{-1} e.
      *   
      * \param i row index
      * \param j column index
      */  
      double p(int i, int j) const; 

      /** 
      * Return sum of elements of the inverse chi matrix.
      */  
      double sumChiInverse() const;

      /**
      * Get number of monomer types.
      */
      int nMonomer() const;

   private:

      // Symmetric matrix of interaction parameters (local copy).
      DMatrix<double> chi_;

      // Inverse of matrix chi_.
      DMatrix<double> chiInverse_;

      // Idempotent matrix P used in residual definition.
      DMatrix<double> p_;

      // Sum of elements of matrix chiInverse_
      double sumChiInverse_;

      /// Number of monomers
      int nMonomer_;

      /// Has private memory been allocated?
      bool isAllocated_;

   };

   // Inline function

   inline int AmbdInteraction::nMonomer() const
   {  return nMonomer_; }

   inline double AmbdInteraction::chi(int i, int j) const
   {  return chi_(i, j); }

   inline double AmbdInteraction::chiInverse(int i, int j) const
   {  return chiInverse_(i, j); }

   inline double AmbdInteraction::p(int i, int j) const
   {  return p_(i, j); }

   inline double AmbdInteraction::sumChiInverse() const
   {  return sumChiInverse_; }

} // namespace Pscf
#endif
