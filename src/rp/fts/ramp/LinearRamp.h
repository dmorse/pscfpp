#ifndef RP_LINEAR_RAMP_H
#define RP_LINEAR_RAMP_H

#include <rp/fts/ramp/Ramp.h>         // base class template
#include <util/containers/DArray.h>   // member
#include <iostream>

// Forward declaration
namespace Pscf {
   namespace Rp {
      template <int D, class T> class RampParameter;
   }
}
namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Linear ramp - parameters vary linearly with step index.
   * 
   * Specializations of this class template are used as base classes for 
   * two closely analogous class templates, both also named LinearRamp,
   * that are defined in Rpc and Rpg namespaces for use in the pscf_rpc
   * and pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *   - D : dimension
   *   - T : Types class, Rpc::Types<D> or Rpg::Types<D>
   *
   * \see \ref rp_LinearRamp_page "Manual Page (Linear Ramp)"
   * \see \ref psfts_ramp_page "Manual Page" (Ramp)"
   * \ingroup Rp_Fts_Ramp_Module
   */
   template <int D, class T>
   class LinearRamp : public Ramp<D,T>
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      */
      LinearRamp(Simulator<D,T>& simulator);

      /**
      * Destructor.
      */
      ~LinearRamp() = default;

      /**
      * Read parameters from parameter file input stream.
      *
      * \param in input parameter file stream
      */
      void readParameters(std::istream& in) override;

      /**
      * Set nStep and complete initialization.
      *
      * This method is called just before the beginning of the main
      * simulation loop.
      *
      * \param nStep number of steps planned for this simulation
      */
      void setup(int nStep) override;

      /**
      * Set new parameters values in associated System and Simulator.
      * 
      * \param iStep  current simulation step index
      */
      void setParameters(int iStep) override;
      
      /**
      * Output information at the end of a simulation.
      */
      void output() override;

   protected:

      using RampT = Ramp<D,T>;
      using RampParameterT = RampParameter<D,T>;

   private:

      // Number of variable parameters
      int nParameter_;

      // Array of variable parameters
      DArray< RampParameter<D,T> > parameters_;

   };

}
}
#endif
