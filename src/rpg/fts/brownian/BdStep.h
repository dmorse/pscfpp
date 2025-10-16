#ifndef RPG_BD_STEP_H
#define RPG_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>
//#include <pscf/cuda/CudaRandom.h>

namespace Pscf {
   class CudaRandom;
   namespace Rpg {
      template <int D> class System;
      template <int D> class BdSimulator;
   }
}

namespace Pscf {
namespace Rpg {

   using namespace Util;


   /**
   * BdStep is an abstract base class for Brownian dynamics steps.
   *
   * The virtual step() method must generate a single step.
   *
   * \ingroup Rpg_Fts_Brownian_Module
   */
   template <int D>
   class BdStep : public ParamComposite
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent BdSimulator<D> object
      */
      BdStep(BdSimulator<D>& simulator);

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      virtual ~BdStep();

      /**
      * Read required parameters from file.
      *
      * Empty default implementation.
      */
      virtual void readParameters(std::istream &in);

      /**
      * Setup before the beginning of each simulation run.
      */
      virtual void setup();

      /**
      * Take a single Brownian dynamics step.
      * 
      * \return true if converged, false if failed to converge.
      */
      virtual bool step() = 0;
      
      /**
      * Decide whether cc fields need to be saved for move.
      * The default implementation is false.
      */
      virtual bool needsCc()
      {  return false; }
      
      /**
      * Decide whether dc fields need to be saved for move.
      * The default implementation is false.
      */
      virtual bool needsDc()
      { return true; }
      
      /**
      * Log output timing results 
      */
      virtual void outputTimers(std::ostream& out);
      
      /**
      * Clear timers 
      */
      virtual void clearTimers();

      // Accessor Functions

      /**
      * Output statistics for this move (at the end of simulation)
      */
      virtual void output();

   protected:

      /**
      * Get parent System object.
      */
      System<D>& system();

      /**
      * Get parent BdSimulator object.
      */
      BdSimulator<D>& simulator();

      /**
      * Get Random number generator of parent System.
      */
      CudaRandom& cudaRandom();

   private:

      /// Pointer to parent BdSimulator object
      BdSimulator<D>* simulatorPtr_;

      /// Pointer to parent System object
      System<D>* systemPtr_;

      /// Pointer to random number generator
      CudaRandom  *cudaRandomPtr_;

   };

   // Protected inline methods

   /*
   * Get parent System object.
   */
   template <int D>
   inline System<D>& BdStep<D>::system()
   {  return *systemPtr_; }

   /*
   * Get parent BdSimulator object.
   */
   template <int D>
   inline BdSimulator<D>& BdStep<D>::simulator()
   {  return *simulatorPtr_; }

   /*
   * Get Random number generator.
   */
   template <int D>
   inline CudaRandom& BdStep<D>::cudaRandom()
   {  return *cudaRandomPtr_; }

   // Explicit instantiation declarations
   extern template class BdStep<1>;
   extern template class BdStep<2>;
   extern template class BdStep<3>;

}
}
#endif
