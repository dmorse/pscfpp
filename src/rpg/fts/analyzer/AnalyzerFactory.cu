/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AnalyzerFactory.h"  

// Subclasses of Analyzer 
#include "TrajectoryWriter.h"
#include "ConcentrationWriter.h"
#include "HamiltonianAnalyzer.h"
#include "BinaryStructureFactor.h"
#include "StepLogger.h"
#include "PerturbationDerivative.h"
#include "BinaryChiDerivative.h"
#include "ConcentrationDerivative.h"
#include "MaxOrderParameter.h"
#include "FourthOrderParameter.h"
#include "CubicLengthDerivative.h"

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   AnalyzerFactory<D>::AnalyzerFactory(Simulator<D>& simulator, Rp::System<D, Rpg::Types<D> >& system)
    : simPtr_(&simulator),
      sysPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Analyzer subclass className.
   */
   template <int D>
   Analyzer<D>* AnalyzerFactory<D>::factory(const std::string &className) const
   {
      Analyzer<D>* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      // Try to match classname
      if (className == "TrajectoryWriter") {
         ptr = new TrajectoryWriter<D>(*simPtr_, *sysPtr_);
      } else if (className == "ConcentrationWriter") {
         ptr = 
           new ConcentrationWriter<D>(*simPtr_, *sysPtr_);
      } else if (className == "HamiltonianAnalyzer") {
         ptr = new HamiltonianAnalyzer<D>(*simPtr_, *sysPtr_);
      } else if (className == "BinaryStructureFactor") {
         ptr = new BinaryStructureFactor<D>(*simPtr_, *sysPtr_);
      } else if (className == "StepLogger") {
         ptr = new StepLogger<D>(*simPtr_, *sysPtr_);
      } else if (className == "PerturbationDerivative") {
         ptr = new PerturbationDerivative<D>(*simPtr_, 
                                             *sysPtr_);
      } else if (className == "BinaryChiDerivative") {
         ptr = new BinaryChiDerivative<D>(*simPtr_, *sysPtr_);
      } else if (className == "ConcentrationDerivative") {
         ptr = new ConcentrationDerivative<D>(*simPtr_, 
                                              *sysPtr_);
      } else if (className == "MaxOrderParameter") {
         ptr = new MaxOrderParameter<D>(*simPtr_, *sysPtr_);
      } else if (className == "FourthOrderParameter") {
         ptr = new FourthOrderParameter<D>(*simPtr_, *sysPtr_);
      } else if (className == "CubicLengthDerivative") {
         ptr = new CubicLengthDerivative<D>(*simPtr_, *sysPtr_);
      }

      return ptr;
   }

   // Explicit instantiation definitions
   template class AnalyzerFactory<1>;
   template class AnalyzerFactory<2>;
   template class AnalyzerFactory<3>;

}
}
