/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#if 0
#include "BinaryChiDerivative.h"                    // header
#include <rp/system/System.h>
#include <rp/fts/simulator/Simulator.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/WFields.h>
#include <rp/field/CFields.h>
#endif

#include <rp/fts/analyzer/BinaryChiDerivative.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class BinaryChiDerivative<1,CPT>;
      template class BinaryChiDerivative<2,CPT>;
      template class BinaryChiDerivative<3,CPT>;
   }
}
