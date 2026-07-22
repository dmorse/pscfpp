#ifndef RPC_SWEEP_FACTORY_H
#define RPC_SWEEP_FACTORY_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/sweep/SweepFactory.h>  
#include <pscf/backends/CPT.h>  

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class SweepFactory<1,CPT>;
      extern template class SweepFactory<2,CPT>;
      extern template class SweepFactory<3,CPT>;
   }
}
#endif
