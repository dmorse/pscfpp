/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backends/CPT.h>
#include <rp/fts/trajectory/TrajectoryReaderFactory.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class TrajectoryReaderFactory<1,CPT>;
      template class TrajectoryReaderFactory<2,CPT>;
      template class TrajectoryReaderFactory<3,CPT>;
   }
}
