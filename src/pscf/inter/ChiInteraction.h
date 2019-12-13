#ifndef PSCF_CHI_INTERACTION_H
#define PSCF_CHI_INTERACTION_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Interaction.h"             // base class
#include <util/global.h>                  

namespace Pscf {

   using namespace Util;

   /**
   * Flory-Huggins excess free energy model.
   *
   * \ingroup Pscf_Module
   */
   class ChiInteraction : public  Interaction
   {

   public:

      /**
      * Constructor.
      */
      ChiInteraction();

      /**
      * Destructor.
      */
      virtual ~ChiInteraction();

      /**
      * Read chi parameters.
      *
      * Must be called after setNMonomer.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Compute excess Helmholtz free energy per monomer.
      *
      * \param c array of concentrations, for each type (input)
      */
      virtual 
      double fHelmholtz(Array<double> const & c) const;

      /**
      * Compute chemical potential from concentration.
      *
      * \param c array of concentrations, for each type (input)
      * \param w array of chemical potentials for types (ouptut) 
      */
      virtual 
      void computeW(Array<double> const & c, Array<double>& w) 
      const;

      /**
      * Compute concentration from chemical potential fields.
      *
      * \param w array of chemical potentials for types (inut) 
      * \param c array of vol. fractions, for each type (output)
      * \param xi Langrange multiplier pressure (output)
      */
      virtual 
      void computeC(Array<double> const & w, Array<double>& c, double& xi)
      const;

      /**
      * Compute Langrange multiplier xi from chemical potential fields.
      *
      * \param w array of chemical potentials for types (inut) 
      * \param xi Langrange multiplier pressure (output)
      */
      virtual 
      void computeXi(Array<double> const & w, double& xi)
      const;

      /**
      * Compute second derivatives of free energy.
      *
      * Upon return, the elements of the square matrix dWdC, are
      * given by derivatives dWdC(i,j) = dW(i)/dC(j), which are
      * also second derivatives of the interaction free energy. 
      * For this Flory-Huggins chi parameter model, this is simply 
      * given by the chi matrix dWdC(i,j) = chi(i, j).
      *
      * \param c array of concentrations, for each type (input)
      * \param dWdC matrix of derivatives (output) 
      */
      virtual 
      void computeDwDc(Array<double> const & c, Matrix<double>& dWdC)
      const;

      /**
      * Return one element of the chi matrix.
      *
      * \param i row index
      * \param j column index
      */
      double chi(int i, int j);

      /**
      * Return one element of the inverse chi matrix.
      *
      * \param i row index
      * \param j column index
      */
      double chiInverse(int i, int j);

      /** 
      * Return one element of the idempotent matrix.
      *   
      * \param i row index
      * \param j column index
      */  
      double idemp(int i, int j); 

      double sum_inv();
   private:

      DMatrix<double> chi_;

      DMatrix<double> chiInverse_;

      DMatrix<double> idemp_;

      double sum_inv_;

   };

   // Inline function

   inline double ChiInteraction::chi(int i, int j)
   {  return chi_(i, j); }

   inline double ChiInteraction::chiInverse(int i, int j)
   {  return chiInverse_(i, j); }

   inline double ChiInteraction::idemp(int i, int j)
   {  return idemp_(i, j); }

   inline double ChiInteraction::sum_inv()
   {  return sum_inv_; }

} // namespace Pscf
#endif
