#ifndef RPC_INTRACORRELATION_H
#define RPC_INTRACORRELATION_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/compressor/IntraCorrelation.h>  // class template
#include <pscf/cpu/CppTp.h>                    // class argument
#include <prdc/field/cpu/FftwDRArray.h>          // class member

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class IntraCorrelation<1, CppTp<1> >;
      extern template class IntraCorrelation<2, CppTp<2> >;
      extern template class IntraCorrelation<3, CppTp<3> >;
   }
}
#endif
