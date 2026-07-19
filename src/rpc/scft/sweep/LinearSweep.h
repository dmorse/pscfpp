#ifndef RPC_LINEAR_SWEEP_H
#define RPC_LINEAR_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/sweep/LinearSweep.h>      // base class template
#include <pscf/cpu/Cpp.h>               // base class argument
#include <rpc/scft/sweep/Sweep.h>           // indirect base class
#include <rpc/scft/sweep/SweepParameter.h>  // indirect base member


// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class LinearSweep<1, Cpp<1> >;
      extern template class LinearSweep<2, Cpp<2> >;
      extern template class LinearSweep<3, Cpp<3> >;
   }
}
#endif
