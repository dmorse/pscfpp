/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AnalyzerFactory.h"

// Subclasses of Analyzer
#include <rp/fts/analyzer/StepLogger.h>
#include <rp/fts/analyzer/TrajectoryWriter.h>
#include <rp/fts/analyzer/ConcentrationWriter.h>
#include <rp/fts/analyzer/HamiltonianAnalyzer.h>
#include <rp/fts/analyzer/BinaryStructureFactor.h>
#include <rp/fts/analyzer/MaxOrderParameter.h>
#include <rp/fts/analyzer/FourthOrderParameter.h>
#include <rp/fts/analyzer/BinaryChiDerivative.h>
#include <rp/fts/analyzer/CubicLengthDerivative.h>
#include <rp/fts/analyzer/ConcentrationDerivative.h>
#include <rp/fts/analyzer/PerturbationDerivative.h>

#include <rp/fts/analyzer/AnalyzerFactory.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class AnalyzerFactory<1,CPT>;
      template class AnalyzerFactory<2,CPT>;
      template class AnalyzerFactory<3,CPT>;
   }
}
