/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/backends/CPT.h>
#include <rp/fts/analyzer/TrajectoryWriter.tpp> 

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class TrajectoryWriter<1,CPT>;
      template class TrajectoryWriter<2,CPT>;
      template class TrajectoryWriter<3,CPT>;
   }
}

