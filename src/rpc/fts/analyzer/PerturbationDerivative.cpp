/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PerturbationDerivative.h"
#include <rpc/fts/perturbation/Perturbation.h>
#include <rp/fts/analyzer/PerturbationDerivative.tpp>
#include <rpc/system/System.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>
#include <rpc/fts/simulator/Simulator.h>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class PerturbationDerivative< 1, Cpp<1> >;
      template class PerturbationDerivative< 2, Cpp<2> >;
      template class PerturbationDerivative< 3, Cpp<3> >;
   }
}
