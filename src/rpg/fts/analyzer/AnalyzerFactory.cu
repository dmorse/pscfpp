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
   AnalyzerFactory<D>::AnalyzerFactory(
                         Rp::Simulator<D, Rpg::Types<D> >& simulator, 
                         Rp::System<D, Rpg::Types<D> >& system)
    : simPtr_(&simulator),
      sysPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Analyzer subclass className.
   */
   template <int D>
   Rp::Analyzer<D, Rpg::Types<D> >* 
   AnalyzerFactory<D>::factory(const std::string &className) const
   {
      Rp::Analyzer<D, Rpg::Types<D> >* ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      // Try to match classname
      if (className == "TrajectoryWriter") {
         ptr = new Rp::TrajectoryWriter<D, Rpg::Types<D> >(*simPtr_, *sysPtr_);
      } else if (className == "ConcentrationWriter") {
         ptr = new Rp::ConcentrationWriter<D, Rpg::Types<D> >(*simPtr_, *sysPtr_);
      } else if (className == "HamiltonianAnalyzer") {
         ptr = new Rp::HamiltonianAnalyzer<D, Rpg::Types<D> >(*simPtr_, *sysPtr_);
      } else if (className == "BinaryStructureFactor") {
         ptr = new BinaryStructureFactor<D>(*simPtr_, *sysPtr_);
      } else if (className == "StepLogger") {
         ptr = new Rp::StepLogger<D, Rpg::Types<D> >(*simPtr_, *sysPtr_);
      } else if (className == "PerturbationDerivative") {
         ptr = new Rp::PerturbationDerivative<D, Rpg::Types<D> >(*simPtr_, *sysPtr_);
      } else if (className == "BinaryChiDerivative") {
         ptr = new Rp::BinaryChiDerivative<D, Rpg::Types<D> >(*simPtr_, *sysPtr_);
      } else if (className == "ConcentrationDerivative") {
         ptr = new Rp::ConcentrationDerivative<D, Rpg::Types<D> >(*simPtr_, *sysPtr_);
      } else if (className == "MaxOrderParameter") {
         ptr = new MaxOrderParameter<D>(*simPtr_, *sysPtr_);
      } else if (className == "FourthOrderParameter") {
         ptr = new FourthOrderParameter<D>(*simPtr_, *sysPtr_);
      } else if (className == "CubicLengthDerivative") {
         ptr = new Rp::CubicLengthDerivative<D, Rpg::Types<D> >(*simPtr_, *sysPtr_);
      }

      return ptr;
   }

   // Explicit instantiation definitions
   template class AnalyzerFactory<1>;
   template class AnalyzerFactory<2>;
   template class AnalyzerFactory<3>;

}
}
