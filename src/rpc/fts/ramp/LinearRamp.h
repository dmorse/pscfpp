#ifndef RPC_LINEAR_RAMP_H
#define RPC_LINEAR_RAMP_H

#include <rpc/fts/ramp/Ramp.h>           // base class
#include <rpc/fts/ramp/RampParameter.h>  // member (templ parameter)
#include <util/containers/DArray.h>      // member (template)

namespace Pscf {
namespace Rpc {

   using namespace Util;

   template <int D> class Simulator;

   /**
   * Linear ramp - parameters vary linearly with step index.
   * 
   * \see
   * <ul>
   *   <li> \ref rpc_LinearRamp_page "Manual Page" </li>
   *   <li> \ref psfts_ramp_page "Manual Page" </li>
   *   <li> Ramp </li>
   * </ul>
   *
   * \ingroup Rpc_Fts_Ramp_Module
   */
   template <int D>
   class LinearRamp : public Ramp<D>
   {

   public:

      /**
      * Constructor.
      */
      LinearRamp(Simulator<D>& simulator);

      /**
      * Destructor.
      */
      virtual ~LinearRamp();

      /**
      * Read parameters from parameter file input stream.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Set nStep and complete initialization.
      *
      * This method is called just before the beginning of the main
      * simulation loop.
      *
      * \param nStep number of steps planned for this simulation
      */
      virtual void setup(int nStep);

      /**
      * Set new parameters values in associated System and Simulator.
      * 
      * \param iStep  current simulation step index
      */
      virtual void setParameters(int iStep);
      
      /**
      * Output initial and final parameter values of linear ramp
      * at the end of the simulation.
      */
      virtual void output();

   protected:

      using Ramp<D>::nStep_;
      using Ramp<D>::simulator;

   private:

      // Number of variable parameters
      int nParameter_;

      // Array of variable parameters
      DArray< RampParameter<D> > parameters_;

   };


   // Explicit instantiation declarations
   extern template class LinearRamp<1>;
   extern template class LinearRamp<2>;
   extern template class LinearRamp<3>;

}
}
#endif
