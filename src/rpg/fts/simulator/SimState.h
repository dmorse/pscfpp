#ifndef RPG_SIM_STATE_H
#define RPG_SIM_STATE_H

/*
*PSCF - Polymer Self-Consistent Field
*
*Copyright 2015 - 2025, The Regents of the University of Minnesota
*Distributed under the terms of the GNU General Public License.
*/

#include<rp/fts/simulator/SimState.h>   // base class template
#include<prdc/cuda/RField.h>            // base class template argument

namespace Pscf {
namespace Rpg {

   using namespace Util;
   
   /**
   * SimState stores the state used by an FTS simulation.
   *
   * This class is used to restore the state of FTS simulation after 
   * an attempted move or step that is rejected or fails to converge. 
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of base class template Rp::SimState, 
   * and inherit their public interface and almost all of their source
   * code from this base class.  
   *
   * \see Rp::SimState
   * \ingroup Rpg_Fts_Simulator_Module
   */
   template <int D>
   struct SimState : Rp::SimState< D, Prdc::Cuda::RField<D> >
   {};

}
}

//Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template struct SimState<1, Prdc::Cuda::RField<1> >;
      extern template struct SimState<2, Prdc::Cuda::RField<2> >;
      extern template struct SimState<3, Prdc::Cuda::RField<3> >;
   }
   namespace Rpg {
      extern template struct SimState<1>;
      extern template struct SimState<2>;
      extern template struct SimState<3>;
   }
}
#endif
