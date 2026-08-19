/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backends/CUT.h>
#include <rp/fts/trajectory/TrajectoryReaderFactory.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class TrajectoryReaderFactory<1,CUT>;
      template class TrajectoryReaderFactory<2,CUT>;
      template class TrajectoryReaderFactory<3,CUT>;
   }
}
