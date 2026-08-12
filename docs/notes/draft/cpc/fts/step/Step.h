#ifndef CPC_BD_STEP_H
#define CPC_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include <util/param/ParamComposite.h>
#include <util/random/Random.h>
#include <util/global.h>

namespace Pscf {
namespace Cpc {

   using namespace Util;

   template <int D> class System;
   template <int D> class Simulator;

   /**
   * Step is an abstract base class for Brownian dynamics steps.
   *
   * The virtual step() method must generate a single step.
   *
   * \ingroup Cpc_Fts_Step_Module
   */
   template <int D>
   class Step : public ParamComposite
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent Simulator<D> object
      */
      Step(Simulator<D>& simulator);

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      virtual ~Step();

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
      * Get parent Simulator object.
      */
      Simulator<D>& simulator();

      /**
      * Get Random number generator of parent Simulator.
      */
      Random& random();

   private:
      
      /// Pointer to parent Simulator object
      Simulator<D>* simulatorPtr_;

      /// Pointer to parent System object
      System<D>* systemPtr_;

      /// Pointer to random number generator of parent Simulator
      Random  *randomPtr_;

   };

   // Protected inline methods
   
   /*
   * Get parent System object.
   */
   template <int D>
   inline System<D>& Step<D>::system()
   {  return *systemPtr_; }

   /*
   * Get parent Simulator object.
   */
   template <int D>
   inline Simulator<D>& Step<D>::simulator()
   {  return *simulatorPtr_; }

   /*
   * Get Random number generator.
   */
   template <int D>
   inline Random& Step<D>::random()
   {  return *randomPtr_; }

   // Explicit instantiation declarations
   extern template class Step<1>;
   extern template class Step<2>;
   extern template class Step<3>;

}
}
#endif
