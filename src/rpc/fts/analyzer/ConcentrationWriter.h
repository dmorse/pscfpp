#ifndef RPC_CONCENTRATION_WRITER_H
#define RPC_CONCENTRATION_WRITER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/ConcentrationWriter.h>  // base class template
#include <pscf/cpu/CppTp.h>                     // base class argument
#include <rpc/fts/analyzer/Analyzer.h>            // indirect base class

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      extern template class ConcentrationWriter<1, CppTp<1> >;
      extern template class ConcentrationWriter<2, CppTp<2> >;
      extern template class ConcentrationWriter<3, CppTp<3> >;
   }
}
#endif
