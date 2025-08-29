#ifndef RPC_RAMP_H
#define RPC_RAMP_H

#include <util/param/ParamComposite.h>      // base class

namespace Pscf {
namespace Rpc {

   using namespace Util;

   template <int D> class Simulator;

   /**
   * Class that varies parameters during a simulation (abstract).
   *
   * \see \ref psfts_ramp_page "Manual Page"
   *
   * \ingroup Rpc_Fts_Ramp_Module
   */
   template <int D>
   class Ramp : public ParamComposite
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent Simulator
      */
      Ramp(Simulator<D>& simulator);

      /**
      * Destructor.
      */
      virtual ~Ramp();

      /**
      * Final setup before simulation loop, set value of nStep.
      *
      * This method must be called just before the beginning of the main
      * simulation loop, after an initial configuration is known. It must
      * set and store the value of nStep (the number of steps planned for
      * the simulation) and complete any initialization that cannot be
      * completed in the readParam method.
      *
      * The default implementation simply stores a value of nStep.
      *
      * \param nStep number of steps planned for this simulation
      */
      virtual void setup(int nStep);

      /**
      * Set new parameters values in associated System and Simulator.
      * 
      * \param iStep  current simulation step index
      */
      virtual void setParameters(int iStep) = 0;
      
      /**
      * Output any results at the end of the simulation.
      *
      * The default implementation is an empty function.
      */
      virtual void output()
      {}

      /**
      * Get parent Simulator<D> by const reference.
      */
      Simulator<D> const & simulator() const;

   protected:

      /**
      * Get parent Simulator<D> by non-const reference.
      */
      Simulator<D>& simulator();

      /// Number of steps planned for this simulation (set in setup).
      int nStep_;

   private:

      /// Pointer to parent Simulator.
      Simulator<D>* simulatorPtr_;

   };

   // Inline methods

   // Return parent simulator by const reference.
   template <int D>
   inline Simulator<D> const & Ramp<D>::simulator() const
   {
      assert(simulatorPtr_);  
      return *simulatorPtr_; 
   }

   // Return parent simulator by non-const reference.
   template <int D>
   inline Simulator<D> & Ramp<D>::simulator() 
   {  
      assert(simulatorPtr_);
      return *simulatorPtr_; 
   }

   // Explicit instantiation declarations
   extern template class Ramp<1>;
   extern template class Ramp<2>;
   extern template class Ramp<3>;

}
}
#endif
