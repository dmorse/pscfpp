/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Solvent.h"

namespace Pscf { 
namespace Cyln
{ 

   /*
   * Constructor
   */
   Solvent::Solvent()
   {  setClassName("Solvent"); }

   /*
   * Destructor
   */
   Solvent::~Solvent()
   {}

   /*
   * Compute monomer concentration field and mu or phi.
   */
   void Solvent::compute(WField const & wField)
   {
   }

}
}
