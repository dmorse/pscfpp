#ifndef PSCF_FH_INTERACTION_H
#define PSCF_FH_INTERACTION_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/containers/DMatrix.h>      // member (template)
#include <util/containers/Matrix.h>       // argument (template)
#include <util/containers/Array.h>        // argument (template)
#include <util/global.h>                  

// Forward declaration
namespace Pscf {
   class Interaction;
}

namespace Pscf {

   using namespace Util;

   /**
   * Flory-Huggins interaction model.
   *
   * \ingroup Pscf_FloryHuggins_Module
   */
   class FhInteraction : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      FhInteraction();

      /**
      * Copy constructor. 
      */
      FhInteraction(FhInteraction const & other);

      /**
      * Constructor - copy from Interaction.
      */
      FhInteraction(Interaction const & other);

      /**
      * Destructor.
      */
      virtual ~FhInteraction();

      /**
      * Assignment.
      * 
      * \param other  rhs of assignment
      */
      FhInteraction& operator = (FhInteraction const & other);

      /**
      * Assignment from Interaction.
      *
      * \param other  rhs of assignment
      */
      FhInteraction& operator = (Interaction const & other);

      /**
      * Set the number of monomer types.
      *
      * \param nMonomer number of monomer types.
      */
      void setNMonomer(int nMonomer);

      /**
      * Read chi parameters.
      *
      * Must be called after setNMonomer.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Assign values for all elements of the chi matrix.
      *
      * On entry:
      *
      *   - nMonomer must be set to a positive value
      *   - The input chi matrix must be nMonomer x nMonomer
      *   - The input chi matrix must be a symmetric 
      *
      * \param chi  input chi matrix
      */
      void setChi(Matrix<double> const & chi);

      /**
      * Assign one element of the chi matrix.
      *
      * \param i row index
      * \param j column index
      * \param chi  input value of chi
      */
      void setChi(int i, int j, double chi);

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
      * \param w array of chemical potentials for types (input) 
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
      * \param c  array of concentrations, for each type (input)
      * \param dWdC  matrix of derivatives (output) 
      */
      virtual 
      void computeDwDc(Array<double> const & c, Matrix<double>& dWdC)
      const;

      /**
      * Return the chi matrix by const reference.
      */
      Matrix<double> const & chi() const;

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
      * Get number of monomer types.
      */
      int nMonomer() const;

   private:

      // Symmetric matrix of interaction parameters.
      DMatrix<double> chi_;

      // Inverse of matrix chi_.
      DMatrix<double> chiInverse_;

      /// Number of monomers
      int nMonomer_;

      /**
      * Compute the inverse of the chi matrix.
      * Must be called after making any changes to the chi matrix.
      */
      void updateMembers();

      void setChiZero();

   };

   // Inline function

   inline int FhInteraction::nMonomer() const
   {  return nMonomer_; }

   inline Matrix<double> const &  FhInteraction::chi() const
   {  return chi_; }

   inline double FhInteraction::chi(int i, int j) const
   {  return chi_(i, j); }

   inline double FhInteraction::chiInverse(int i, int j) const
   {  return chiInverse_(i, j); }

} // namespace Pscf
#endif
