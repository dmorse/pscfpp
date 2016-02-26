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
      * Compute wField from monomer cFields and pressure.
      */
      virtual void computeWField(Array<double> const & cField, double pressure, 
                                Array<double>& wField);

      /**
      * Compute cField from wField. 
      */
      virtual void computecField(Array<double> const & wField, 
                                 Array<double>& cField);

   private:

      DMatrix<double> chi_;

      DMatrix<double> chiInverse_;

   };

} // namespace Pscf
#endif
