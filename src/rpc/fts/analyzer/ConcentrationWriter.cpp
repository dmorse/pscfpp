/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ConcentrationWriter.h"

#include <rp/fts/simulator/Simulator.h>
#include <rp/system/System.h>
#include <rp/solvers/Mixture.h>
#include <rp/field/Domain.h>
#include <rp/field/FieldIo.h>
#include <rp/field/CFields.h>
#include <rp/field/WFields.h>
#include <prdc/field/cpu/RField.h>

#include <rp/fts/analyzer/ConcentrationWriter.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class ConcentrationWriter<1,CPT>;
      template class ConcentrationWriter<2,CPT>;
      template class ConcentrationWriter<3,CPT>;
   }
}
