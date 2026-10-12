#ifndef RP_LINEAR_RAMP_H
#define RP_LINEAR_RAMP_H

#include <rp/fts/ramp/Ramp.h>          // base class template
#include <rp/fts/ramp/RampParameter.h> // member
#include <util/containers/DArray.h>    // member

#include <pscf/backend/TmplDeclare.h>
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
   * Linear ramp - on or more parameters vary linearly with step index.
   * 
   * Template parameters:
   *
   *   - D : dimension
   *   - T : backend identifier class (CPT or CUT)
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

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(LinearRamp)

} // namespace Rp
} // namespace Pscf
#endif
