/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Perturbation.h"
#include <rpc/fts/simulator/Simulator.h>
#include <prdc/field/cpu/RField.h>

#include <rp/fts/perturbation/Perturbation.tpp> // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class Perturbation<1, Cpp<1> >;
      template class Perturbation<2, Cpp<2> >;
      template class Perturbation<3, Cpp<3> >;
   }
}
