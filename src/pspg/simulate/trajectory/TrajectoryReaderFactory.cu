/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryReaderFactory.tpp"

namespace Pscf {
namespace Pspg {
   template class TrajectoryReaderFactory<1>;
   template class TrajectoryReaderFactory<2>;
   template class TrajectoryReaderFactory<3>;
}
}
