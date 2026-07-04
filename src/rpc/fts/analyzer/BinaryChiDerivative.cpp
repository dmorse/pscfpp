/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BinaryChiDerivative.h"                    // header
#include <rpc/system/System.h>
#include <rpc/fts/simulator/Simulator.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/WFields.h>
#include <rpc/field/CFields.h>

#include <rp/fts/analyzer/BinaryChiDerivative.tpp>  // base class implementation

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BinaryChiDerivative< 1, Rpc::Types<1> >;
      template class BinaryChiDerivative< 2, Rpc::Types<2> >;
      template class BinaryChiDerivative< 3, Rpc::Types<3> >;
   }
}
