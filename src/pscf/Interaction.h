#ifndef PSCF_INTERACTION_H
#define PSCF_INTERACTION_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2013, David Morse (morse012@.umn.edu)
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
      * Compute omega from and pressure p.
      *
      * \return error code: 0 for success, 1 for failure.
      */
      virtual void computeOmega(Array<double> const & rho, double p, Array<double>& omega) = 0;

      /**
      * Compute omega from and pressure p.
      *
      * \return error code: 0 for success, 1 for failure.
      */
      virtual void computeConcentration(Array<double> const & omega, Array<double>& omega) = 0;

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
