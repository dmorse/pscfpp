/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AnalyzerFactory.h"

// Subclasses of Analyzer
#include <rpc/fts/analyzer/StepLogger.h>
#include <rpc/fts/analyzer/TrajectoryWriter.h>
#include <rpc/fts/analyzer/ConcentrationWriter.h>
#include <rpc/fts/analyzer/HamiltonianAnalyzer.h>
#include <rpc/fts/analyzer/BinaryStructureFactor.h>
#include <rpc/fts/analyzer/MaxOrderParameter.h>
#include <rpc/fts/analyzer/FourthOrderParameter.h>
#include <rpc/fts/analyzer/BinaryChiDerivative.h>
#include <rpc/fts/analyzer/CubicLengthDerivative.h>
#include <rpc/fts/analyzer/ConcentrationDerivative.h>
#include <rpc/fts/analyzer/PerturbationDerivative.h>

#include <rp/fts/analyzer/AnalyzerFactory.tpp>

// Explicit instantiation definitions
namespace Pscf {
   namespace Rp {
      template class AnalyzerFactory<1, CppTp<1> >;
      template class AnalyzerFactory<2, CppTp<2> >;
      template class AnalyzerFactory<3, CppTp<3> >;
   }
}
