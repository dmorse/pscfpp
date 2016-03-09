#ifndef PSCF_INTERACTION_H
#define PSCF_INTERACTION_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/containers/Array.h>        // argument
#include <util/global.h>                  

namespace Pscf {

   using namespace Util;

   /**
   * Base class for excess free energy models.
   *
   * \ingroup Pscf_Fd1d_Module
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
      double fHelmholtz(Array<double> const & c)
      {  return 0.0; }

      /**
      * Compute chemical potentials from concentrations and pressure p.
      *
      * \param c array of concentrations, for each type (input)
      * \param p pressure / Lagrange multiplier (input)
      * \param w array of chemical potentials, for each type (output)
      */
      virtual 
      void computeW(Array<double> const & c, double p, Array<double>& w)
      {}

      /**
      * Compute concentration and pressure p from chemical potentials.
      *
      * \param w  array of chemical potentials, for each type (input)
      * \param c  array of concentrations, for each type (output)
      * \param p  pressure / Lagrange multiplier (output)
      */
      virtual 
      void computeC(Array<double> const & w, 
                    Array<double>& c,
                    double p)
      {}

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
