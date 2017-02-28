#ifndef PSSP_SOLVENT_TPP
#define PSSP_SOLVENT_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Solvent.h"

namespace Pscf { 
namespace Pssp { 

   /*
   * Constructor
   */
   template <int D>
   Solvent<D>::Solvent()
   {  setClassName("Solvent"); }

   /*
   * Destructor
   */
   Solvent<D>::~Solvent()
   {}

   /*
   * Compute monomer concentration field and mu or phi.
   */
   void Solvent<D>::compute(WField const & wField)
   {
   }

}
}
#endif
