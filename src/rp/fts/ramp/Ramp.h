#ifndef RP_RAMP_H
#define RP_RAMP_H

#include <util/param/ParamComposite.h>      // base class

namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Class that varies parameters during a simulation (abstract).
   *
   * Specializations of this class template are used as base classes for 
   * two closely analogous class templates, both also named Ramp, that
   * are defined in Rpc and Rpg namespaces and used in the pscf_rpc and
   * pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *   - D : dimension
   *   - T : Types class, Rpc::Types<D> or Rpg::Types<D>
   *
   * \see \ref psfts_ramp_page "Manual Page"
   * \ingroup Rp_Fts_Ramp_Module
   */
   template <int D, class T>
   class Ramp : public ParamComposite
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent Simulator
      */
      Ramp(typename T::Simulator& simulator);

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
      * Get parent typename T::Simulator by const reference.
      */
      typename T::Simulator const & simulator() const;

   protected:

      /**
      * Get parent typename T::Simulator by non-const reference.
      */
      typename T::Simulator& simulator();

      /// Number of steps planned for this simulation (set in setup).
      int nStep_;

   private:

      /// Pointer to parent Simulator.
      typename T::Simulator* simulatorPtr_;

   };

   // Inline methods

   // Return parent simulator by const reference.
   template <int D, class T> inline 
   typename T::Simulator const & Ramp<D,T>::simulator() const
   {
      assert(simulatorPtr_);  
      return *simulatorPtr_; 
   }

   // Return parent simulator by non-const reference.
   template <int D, class T> inline 
   typename T::Simulator & Ramp<D,T>::simulator() 
   {  
      assert(simulatorPtr_);
      return *simulatorPtr_; 
   }

}
}
#endif
