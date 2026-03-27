#ifndef RP_LINEAR_RAMP_H
#define RP_LINEAR_RAMP_H

#include <util/containers/DArray.h>      // member
#include <iostream>

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
   class LinearRamp : public T::Ramp
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator object
      */
      LinearRamp(typename T::Simulator& simulator);

      /**
      * Destructor.
      */
      virtual ~LinearRamp();

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

      using RampT = typename T::Ramp;
      using RampParameterT = typename T::RampParameter;

   private:

      // Number of variable parameters
      int nParameter_;

      // Array of variable parameters
      DArray< RampParameterT > parameters_;

   };

}
}
#endif
