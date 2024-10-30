#ifndef RPG_LINEAR_RAMP_H
#define RPG_LINEAR_RAMP_H

#include <rpg/fts/ramp/Ramp.h>           // base class
#include <rpg/fts/ramp/RampParameter.h>  // member (templ parameter)
#include <util/containers/DArray.h>           // member (template)

namespace Pscf {
namespace Rpg {

   using namespace Util;

   template <int D> class Simulator;

   /**
   * Linear ramp - parameters vary linearly with step index.
   *
   * \ingroup Rpg_Fts_Ramp_Module
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

   private:

      // Number of variable parameters
      int nParameter_;

      // Array of variable parameters
      DArray< RampParameter<D> > parameters_;

   };


   #ifndef RPG_LINEAR_RAMP_TPP
   // Suppress implicit instantiation
   extern template class LinearRamp<1>;
   extern template class LinearRamp<2>;
   extern template class LinearRamp<3>;
   #endif

}
}
#endif
