/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryReader.h"

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      template class TrajectoryReader<1, CppTp<1> >;
      template class TrajectoryReader<2, CppTp<2> >;
      template class TrajectoryReader<3, CppTp<3> >;
   }
}
