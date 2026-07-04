/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AnalyzerFactory.h"

// Subclasses of Analyzer
#include "StepLogger.h"
#include "TrajectoryWriter.h"
#include "ConcentrationWriter.h"
#include "HamiltonianAnalyzer.h"
#include "BinaryStructureFactor.h"
#include "MaxOrderParameter.h"
#include "FourthOrderParameter.h"
#include "BinaryChiDerivative.h"
#include "CubicLengthDerivative.h"
#include "ConcentrationDerivative.h"
#include "PerturbationDerivative.h"

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   AnalyzerFactory<D>::AnalyzerFactory(Rp::Simulator<D, Rpc::Types<D> >& simulator,
                                       Rp::System<D, Rpc::Types<D> >& system)
    : simPtr_(&simulator),
      sysPtr_(&system)
   {}

   /*
   * Return a pointer to a instance of Analyzer subclass className.
   */
   template <int D>
   Rp::Analyzer<D, Rpc::Types<D> >* 
   AnalyzerFactory<D>::factory(const std::string &className) const
   {
      Rp::Analyzer<D, Rpc::Types<D> >* ptr = nullptr;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;


      // Try to match classname
      if (className == "StepLogger") {
         ptr = new Rp::StepLogger<D, Rpc::Types<D> >(*simPtr_, *sysPtr_);
      } else if (className == "TrajectoryWriter") {
         ptr = new Rp::TrajectoryWriter<D, Rpc::Types<D> >(*simPtr_, *sysPtr_);
      } else if (className == "ConcentrationWriter") {
         ptr = new Rp::ConcentrationWriter<D, Rpc::Types<D> >(*simPtr_, *sysPtr_);
      } else if (className == "HamiltonianAnalyzer") {
         ptr = new Rp::HamiltonianAnalyzer<D, Rpc::Types<D> >(*simPtr_, *sysPtr_);
      } else if (className == "BinaryStructureFactor") {
         ptr = new Rp::BinaryStructureFactor<D, Rpc::Types<D> >(*simPtr_, *sysPtr_);
      } else if (className == "PerturbationDerivative") {
         ptr = new Rp::PerturbationDerivative<D, Rpc::Types<D> >(*simPtr_, *sysPtr_);
      } else if (className == "BinaryChiDerivative") {
         ptr = new Rp::BinaryChiDerivative<D, Rpc::Types<D> >(*simPtr_, *sysPtr_);
      } else if (className == "CubicLengthDerivative") {
         ptr = new Rp::CubicLengthDerivative<D, Rpc::Types<D> >(*simPtr_, *sysPtr_);
      } else if (className == "ConcentrationDerivative") {
         ptr = new Rp::ConcentrationDerivative<D, Rpc::Types<D> >(*simPtr_, *sysPtr_);
      } else if (className == "MaxOrderParameter") {
         ptr = new Rp::MaxOrderParameter<D, Rpc::Types<D> >(*simPtr_, *sysPtr_);
      } else if (className == "FourthOrderParameter") {
         ptr = new Rp::FourthOrderParameter<D, Rpc::Types<D> >(*simPtr_, *sysPtr_);
      }

      return ptr;
   }

   // Explicit instantiation definitions
   template class AnalyzerFactory<1>;
   template class AnalyzerFactory<2>;
   template class AnalyzerFactory<3>;

}
}
