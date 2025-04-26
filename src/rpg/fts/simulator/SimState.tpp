#ifndef RPG_SIM_STATE_TPP
#define RPG_SIM_STATE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SimState.h"

namespace Pscf {
namespace Rpg {

   using namespace Util;
   using namespace Pscf::Prdc::Cuda;

   /*
   * Constructor.
   */
   template <int D>
   SimState<D>::SimState() 
     : w(),
       wc(), 
       hamiltonian(0.0),
       idealHamiltonian(0.0),
       fieldHamiltonian(0.0),
       needsCc(false),
       needsDc(false),
       needsHamiltonian(false), 
       hasData(false),
       isAllocated(false)
   {}

   /*
   * Destructor.
   */
   template <int D>
   SimState<D>::~SimState() 
   {}

   /*
   * Allocate memory for w fields.
   */ 
   template <int D>
   void SimState<D>::allocate(int nMonomer, IntVec<D> const & dimensions)
   {
      w.allocate(nMonomer);
      wc.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         w[i].allocate(dimensions);
         wc[i].allocate(dimensions);
      }
      
      if (needsCc){
         cc.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            cc[i].allocate(dimensions);
         }
      }
      if (needsDc){
         dc.allocate(nMonomer-1);
         for (int i = 0; i < nMonomer - 1; ++i) {
            dc[i].allocate(dimensions);
         }
      }
      
      isAllocated = true;
   }

}
}
#endif
