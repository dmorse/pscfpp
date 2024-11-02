#ifndef RPC_LINEAR_RAMP_TPP
#define RPC_LINEAR_RAMP_TPP

#include "LinearRamp.h"
#include <rpc/fts/simulator/Simulator.h>

#include <util/containers/DArray.h>
#include <util/global.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /* 
   * Constructor.
   */
   template <int D>
   LinearRamp<D>::LinearRamp(Simulator<D>& simulator)
    : Ramp<D>(simulator)
   {}

   /* 
   * Destructor.
   */
   template <int D>
   LinearRamp<D>::~LinearRamp()
   {}

   template <int D>
   void LinearRamp<D>::readParameters(std::istream& in)
   {
      // Read the number of ramp parameters, allocate parameters_ array
      ParamComposite::read(in, "nParameter", nParameter_);
      UTIL_CHECK(nParameter_ > 0);
      parameters_.allocate(nParameter_);
      
      // Read array of RampParameters, using << operator for each
      ParamComposite::readDArray(in, "parameters", parameters_, nParameter_);

      // Verify net zero change in volume fractions, if these are swept
      double sum = 0.0;
      for (int i = 0; i < nParameter_; ++i) {
         if (parameters_[i].type() == "phi_polymer" || 
             parameters_[i].type() == "phi_solvent") 
         {
            sum += parameters_[i].change();
         }
      }
      UTIL_CHECK(sum > -0.000001);
      UTIL_CHECK(sum < 0.000001);
   }

   template <int D>
   void LinearRamp<D>::setup(int nStep)
   {
      nStep_ = nStep;

      // Set simulator pointer and initial value for each parameter object
      for (int i = 0; i < nParameter_; ++i) {
         parameters_[i].setSimulator(Ramp<D>::simulator());
         parameters_[i].getInitial();
      }
   }

   template <int D>
   void LinearRamp<D>::setParameters(int iStep)
   {
      // Compute a ramp parameter in the range [0,1]
      double s = double(iStep)/double(nStep_);

      // Update the system parameter values
      double newVal;
      for (int i = 0; i < nParameter_; ++i) {
         newVal = parameters_[i].initial() + s*parameters_[i].change();
         parameters_[i].update(newVal);
      }
   }
   
   template <int D>
   void LinearRamp<D>::output()
   {
      for (int i = 0; i < nParameter_; ++i) {
         
         // Output parameter type
         Log::file()<< "Linear Ramp parameter: " 
                    << parameters_[i].type() << std::endl;
                    
         // Output initial value
         Log::file()<< "The initial value is: "
                    << parameters_[i].initial() << " ";
         
         // Output final value 
         Log::file()<< "The final value is: "
                    << parameters_[i].initial() + parameters_[i].change() 
                    << std::endl;
         
         // Output the magnitude of change           
         Log::file()<< "The change is: "
                    << parameters_[i].change() <<std::endl;
      }
   }

}
}
#endif 
