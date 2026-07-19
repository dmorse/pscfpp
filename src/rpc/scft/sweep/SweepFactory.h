#ifndef RPC_SWEEP_FACTORY_H
#define RPC_SWEEP_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/sweep/SweepFactory.h>  
#include <pscf/cpu/Cpp.h>  

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class SweepFactory<1, Cpp<1> >;
      extern template class SweepFactory<2, Cpp<2> >;
      extern template class SweepFactory<3, Cpp<3> >;
   }
}
#endif
