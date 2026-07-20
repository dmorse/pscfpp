/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ConcentrationWriter.h"

#include <rpc/fts/simulator/Simulator.h>
#include <rpc/system/System.h>
#include <rpc/solvers/Mixture.h>
#include <rpc/field/Domain.h>
#include <rpc/field/FieldIo.h>
#include <rpc/field/CFields.h>
#include <rpc/field/WFields.h>
#include <prdc/field/cpu/RField.h>

#include <rp/fts/analyzer/ConcentrationWriter.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ConcentrationWriter<1, CppTp<1> >;
      template class ConcentrationWriter<2, CppTp<2> >;
      template class ConcentrationWriter<3, CppTp<3> >;
   }
}
