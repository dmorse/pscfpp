/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AnalyzerFactory.h"

// Subclasses of Analyzer
#include <rpg/fts/analyzer/StepLogger.h>
#include <rpg/fts/analyzer/TrajectoryWriter.h>
#include <rpg/fts/analyzer/ConcentrationWriter.h>
#include <rpg/fts/analyzer/HamiltonianAnalyzer.h>
#include <rpg/fts/analyzer/BinaryStructureFactor.h>
#include <rpg/fts/analyzer/MaxOrderParameter.h>
#include <rpg/fts/analyzer/FourthOrderParameter.h>
#include <rpg/fts/analyzer/BinaryChiDerivative.h>
#include <rpg/fts/analyzer/CubicLengthDerivative.h>
#include <rpg/fts/analyzer/ConcentrationDerivative.h>
#include <rpg/fts/analyzer/PerturbationDerivative.h>

#include <rp/fts/analyzer/AnalyzerFactory.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class AnalyzerFactory<1,CUT>;
      template class AnalyzerFactory<2,CUT>;
      template class AnalyzerFactory<3,CUT>;
   }
}
