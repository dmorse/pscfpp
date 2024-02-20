#ifndef RPG_MC_STATE_TPP
#define RPG_MC_STATE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McState.h"

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   McState<D>::McState() 
     : w(),
       wc(), 
       hamiltonian(0.0),
       idealHamiltonian(0.0),
       fieldHamiltonian(0.0), 
       hasData(false),
       isAllocated(false)
   {}

   /*
   * Destructor.
   */
   template <int D>
   McState<D>::~McState() 
   {}

   /*
   * Allocate memory for w fields.
   */ 
   template <int D>
   void McState<D>::allocate(int nMonomer, IntVec<D> const & dimensions)
   {
      w.allocate(nMonomer);
      wc.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w[i].allocate(dimensions);
         wc[i].allocate(dimensions);
      }
      isAllocated = true;
   }

}
}
#endif
