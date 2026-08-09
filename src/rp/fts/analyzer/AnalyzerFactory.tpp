/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/fts/analyzer/AnalyzerFactory.h>

// Subclasses of Analyzer
#include <rp/fts/analyzer/StepLogger.h>
#include <rp/fts/analyzer/TrajectoryWriter.h>
#include <rp/fts/analyzer/ConcentrationWriter.h>
#include <rp/fts/analyzer/HamiltonianAnalyzer.h>
#include <rp/fts/analyzer/BinaryChiDerivative.h>
#include <rp/fts/analyzer/CubicLengthDerivative.h>
#include <rp/fts/analyzer/ConcentrationDerivative.h>
#include <rp/fts/analyzer/PerturbationDerivative.h>
#include <rp/fts/analyzer/BinaryStructureFactor.h>
#include <rp/fts/analyzer/MaxOrderParameter.h>
#include <rp/fts/analyzer/FourthOrderParameter.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D, class T>
   AnalyzerFactory<D,T>::AnalyzerFactory(Simulator<D,T>& simulator,
                                       System<D,T>& system)
    : simPtr_(&simulator),
      sysPtr_(&system)
   {}

   /*
   * Return a pointer to a instance of Analyzer subclass className.
   */
   template <int D, class T>
   Analyzer<D,T>* 
   AnalyzerFactory<D,T>::factory(const std::string &className) const
   {
      Analyzer<D,T>* ptr = nullptr;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;


      // Try to match classname
      if (className == "StepLogger") {
         ptr = new StepLogger<D,T>(*simPtr_, *sysPtr_);
      } else if (className == "TrajectoryWriter") {
         ptr = new TrajectoryWriter<D,T>(*simPtr_, *sysPtr_);
      } else if (className == "ConcentrationWriter") {
         ptr = new ConcentrationWriter<D,T>(*simPtr_, *sysPtr_);
      } else if (className == "HamiltonianAnalyzer") {
         ptr = new HamiltonianAnalyzer<D,T>(*simPtr_, *sysPtr_);
      } else if (className == "BinaryStructureFactor") {
         ptr = new BinaryStructureFactor<D,T>(*simPtr_, *sysPtr_);
      } else if (className == "PerturbationDerivative") {
         ptr = new PerturbationDerivative<D,T>(*simPtr_, *sysPtr_);
      } else if (className == "BinaryChiDerivative") {
         ptr = new BinaryChiDerivative<D,T>(*simPtr_, *sysPtr_);
      } else if (className == "CubicLengthDerivative") {
         ptr = new CubicLengthDerivative<D,T>(*simPtr_, *sysPtr_);
      } else if (className == "ConcentrationDerivative") {
         ptr = new ConcentrationDerivative<D,T>(*simPtr_, *sysPtr_);
      } else if (className == "MaxOrderParameter") {
         ptr = new MaxOrderParameter<D,T>(*simPtr_, *sysPtr_);
      } else if (className == "FourthOrderParameter") {
         ptr = new FourthOrderParameter<D,T>(*simPtr_, *sysPtr_);
      }

      return ptr;
   }

}
}
