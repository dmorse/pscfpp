#ifndef PSCF_CHI_INTERACTION_H
#define PSCF_CHI_INTERACTION_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Interaction.h"             // base class
#include <util/global.h>                  

namespace Pscf {

   using namespace Util;

   /**
   * Flory-Huggins excess free energy model.
   *
   * \ingroup Pscf_Fd1d_Module
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
      * Compute omega from monomer concentrations and pressure.
      */
      virtual void computeOmega(Array<double> const & concentration, double pressure, 
                                Array<double>& omega);

      /**
      * Compute omega from and pressure p.
      *
      * \return error code: 0 for success, 1 for failure.
      */
      virtual void computeConcentration(Array<double> const & omega, Array<double>& omega);

   private:

      DMatrix<double> chi_;

      DMatrix<double> chiInverse_;

   };

} // namespace Pscf
#endif
