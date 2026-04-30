#ifndef RP_LINEAR_RAMP_TPP
#define RP_LINEAR_RAMP_TPP

#include "LinearRamp.h"
#include <util/global.h>

namespace Pscf {
namespace Rp {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D, class T>
   LinearRamp<D,T>::LinearRamp(typename T::Simulator& simulator)
    : RampT(simulator)
   {  ParamComposite::setClassName("LinearRamp"); }

   template <int D, class T>
   void LinearRamp<D,T>::readParameters(std::istream& in)
   {
      // Read the number of ramp parameters, allocate parameters_ array
      ParamComposite::read(in, "nParameter", nParameter_);
      UTIL_CHECK(nParameter_ > 0);
      parameters_.allocate(nParameter_);

      // Read array of RampParameters, using << operator for each
      ParamComposite::readDArray(in, "parameters", parameters_, 
                                 nParameter_);

      // Verify net zero change in volume fractions, if these are swept
      double sum = 0.0;
      for (int i = 0; i < nParameter_; ++i) {
         if (parameters_[i].type() == "phi_polymer" ||
             parameters_[i].type() == "phi_solvent")
         {
            sum += parameters_[i].change();
         }
      }
      UTIL_CHECK(sum > -0.0000001);
      UTIL_CHECK(sum < 0.0000001);
   }

   template <int D, class T>
   void LinearRamp<D,T>::setup(int nStep)
   {
      RampT::nStep_ = nStep;

      // Set simulator pointer and initial value for each parameter object
      for (int i = 0; i < nParameter_; ++i) {
         parameters_[i].setSimulator(RampT::simulator());
         parameters_[i].getInitial();
      }
   }

   template <int D, class T>
   void LinearRamp<D,T>::setParameters(int iStep)
   {
      // Compute a ramp parameter in the range [0,1]
      double s = double(iStep)/double(RampT::nStep_);

      // Update the system parameter values
      double newVal;
      for (int i = 0; i < nParameter_; ++i) {
         newVal = parameters_[i].initial() + s*parameters_[i].change();
         parameters_[i].update(newVal);

         // Update chiEvals and chiEvecs if parameter is chi
         if (parameters_[i].type() == "chi") {
            RampT::simulator().analyzeChi();
         }

      }
   }

   template <int D, class T>
   void LinearRamp<D,T>::output()
   {
      for (int i = 0; i < nParameter_; ++i) {
         Log::file() << "Parameter type: "
                     << parameters_[i].type() << std::endl;
         Log::file() << "Initial: "
                     << parameters_[i].initial() << " ";
         Log::file()<< "Change: "
                    << parameters_[i].change() << std::endl;
      }
   }

}
}
#endif
