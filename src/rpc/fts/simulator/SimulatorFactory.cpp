/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "SimulatorFactory.tpp"

namespace Pscf {
namespace Rpc {
   template class SimulatorFactory<1>;
   template class SimulatorFactory<2>;
   template class SimulatorFactory<3>;
}
}
