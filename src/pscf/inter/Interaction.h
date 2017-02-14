#ifndef PSCF_INTERACTION_H
#define PSCF_INTERACTION_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/containers/Array.h>        // argument (template)
#include <util/containers/Matrix.h>        // argument (template)
#include <util/global.h>                  

namespace Pscf {

   using namespace Util;

   /**
   * Base class for excess free energy models.
   *
   * \ingroup Pscf_Module
   */
   class Interaction : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Interaction();

      /**
      * Destructor.
      */
      virtual ~Interaction();

      /**
      * Set the number of monomer types.
      *
      * \param nMonomer number of monomer types.
      */
      void setNMonomer(int nMonomer);

      /**
      * Compute excess Helmholtz free energy per monomer.
      *
      * \param c array of concentrations, for each type (input)
      */
      virtual 
      double fHelmholtz(Array<double> const & c) const = 0;
      // { return 0.0; }

      /**
      * Compute interaction contributions to chemical potentials.
      *
      * The resulting chemical potential fields are those obtained
      * with a vanishing Lagrange multiplier / pressure field, xi = 0.
      *
      * \param c array of concentrations, for each type (input)
      * \param w array of chemical potentials, for each type (output)
      */
      virtual 
      void computeW(Array<double> const & c, Array<double>& w) 
      const = 0;
      // {}

      /**
      * Compute concentration and xi from chemical potentials.
      *
      * \param w  array of chemical potentials, for each type (input)
      * \param c  array of concentrations, for each type (output)
      * \param xi  Lagrange multiplier pressure (output)
      */
      virtual 
      void computeC(Array<double> const & w, 
                    Array<double>& c, double& xi) const = 0;
      // {}

      /**
      * Compute Langrange multiplier xi from chemical potentials.
      *
      * \param w  array of chemical potentials, for each type (input)
      * \param xi  Lagrange multiplier pressure (output)
      */
      virtual 
      void computeXi(Array<double> const & w, double& xi) const = 0;
      //{}

      /**
      * Compute matrix of derivatives of w fields w/ respect to c fields.
      *
      * Upon return, the elements of the matrix dWdC are given by 
      * derivatives elements dWdC(i,j) = dW(i)/dC(j), which are also
      * second derivatives of fHelmholtz with respect to concentrations.
      *
      * \param c array of concentrations, for each type (input)
      * \param dWdC square symmetric matrix of derivatives (output)
      */
      virtual 
      void computeDwDc(Array<double> const & c, Matrix<double>& dWdC)
      const = 0;
      //{}

      /**
      * Get number of monomer types.
      */
      int nMonomer() const;

   private:

      /// Number of monomers
      int nMonomer_;

   };

   /**
   * Get number of monomer types.
   */
   inline 
   int Interaction::nMonomer() const
   {  return nMonomer_; }

} // namespace Pscf
#endif
