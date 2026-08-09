#ifndef RPC_BD_SIMULATOR_H
#define RPC_BD_SIMULATOR_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/brownian/BdSimulator.h>       // base class template
#include <pscf/backends/CPT.h>                  // template argument
#include <rp/fts/simulator/Simulator.h>       // indirect base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class BdSimulator<1,CPT>;
      extern template class BdSimulator<2,CPT>;
      extern template class BdSimulator<3,CPT>;
   }
}
#endif
