/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryReaderFactory.tpp"

namespace Pscf {
namespace Rpc {
   template class TrajectoryReaderFactory<1>;
   template class TrajectoryReaderFactory<2>;
   template class TrajectoryReaderFactory<3>;
}
}
