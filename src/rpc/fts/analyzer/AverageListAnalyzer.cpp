/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rpc/fts/analyzer/AverageListAnalyzer.h>
#include <rpc/system/System.h>

#include <rp/fts/analyzer/AverageListAnalyzer.tpp> // implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class AverageListAnalyzer<1, Rpc::Types<1> >;
      template class AverageListAnalyzer<2, Rpc::Types<2> >;
      template class AverageListAnalyzer<3, Rpc::Types<3> >;
   }
}
