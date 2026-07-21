#ifndef RPG_SIM_STATE_H
#define RPG_SIM_STATE_H

/*
*PSCF - Polymer Self-Consistent Field
*
*Copyright 2015 - 2025, The Regents of the University of Minnesota
*Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/simulator/SimState.h>  // class template
#include <pscf/cuda/CudaTp.h>           // template argument

//Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class SimState<1, CudaTp<1> >;
      extern template class SimState<2, CudaTp<2> >;
      extern template class SimState<3, CudaTp<3> >;
   }
}
#endif
