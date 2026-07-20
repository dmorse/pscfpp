#ifndef RPC_SWEEP_H
#define RPC_SWEEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/scft/sweep/Sweep.h>             // base class template
#include <pscf/cpu/CppTp.h>                // base class argument
#include <rpc/scft/sweep/BasisFieldState.h>  // indirect base argument

// Explicit instantiation declarations
namespace Pscf {
   extern template class SweepTmpl< Rp::BasisFieldState<1, CppTp<1> > >;
   extern template class SweepTmpl< Rp::BasisFieldState<2, CppTp<2> > >;
   extern template class SweepTmpl< Rp::BasisFieldState<3, CppTp<3> > >;
   namespace Rp {
      extern template class Sweep<1, CppTp<1> >;
      extern template class Sweep<2, CppTp<2> >;
      extern template class Sweep<3, CppTp<3> >;
   }
}
#endif
