#ifndef RPC_BD_STEP_H
#define RPC_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include <util/param/ParamComposite.h>
#include <util/random/Random.h>
#include <util/global.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   template <int D> class System;
   template <int D> class BdSimulator;

   /**
   * BdStep is an abstract base class for Brownian dynamics steps.
   *
   * The virtual step() method must generate a single step.
   *
   * \ingroup Rpc_Fts_Brownian_Module
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
      *
      * \param in input stream from which to read
      */
      virtual void readParameters(std::istream &in);

      /**
      * Setup before the beginning of each simulation run.
      */
      virtual void setup();

      /**
      * Take a single Brownian dynamics step.
      * 
      * \return true if the compressor converged, false if it failed.
      */
      virtual bool step() = 0;
      
      /**
      * Do cc concentration components need to be saved before a step?
      *
      * The default implementation returns false.
      *
      * \return true to save, or false otherwise
      */
      virtual bool needsCc()
      {  return false; }
      
      /**
      * Do dc derivative components need to be saved before a step?
      *
      * The default implementation returns false.
      *
      * \return true to save, or false otherwise
      */
      virtual bool needsDc()
      { return true; }
      
      /**
      * Output timing results to ostream.
      *
      * \param out output stream
      */
      virtual void outputTimers(std::ostream& out);
      
      /**
      * Clear timers. 
      */
      virtual void clearTimers();
      
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
      Random& random();

   private:
      
      /// Pointer to parent BdSimulator object
      BdSimulator<D>* simulatorPtr_;

      /// Pointer to parent System object
      System<D>* systemPtr_;

      /// Pointer to random number generator
      Random  *randomPtr_;

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
   inline Random& BdStep<D>::random()
   {  return *randomPtr_; }

   #ifndef RPC_BD_STEP_TPP
   // Suppress implicit instantiation
   extern template class BdStep<1>;
   extern template class BdStep<2>;
   extern template class BdStep<3>;
   #endif

}
}
#endif
