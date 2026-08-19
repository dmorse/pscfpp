#ifndef RP_BD_STEP_H
#define RP_BD_STEP_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // base class
#include <util/global.h>
#include <pscf/backends/TmplDeclare.h>

// Forward declaration
namespace Util {
   class Random;
}
namespace Pscf {
   namespace Rp {
      template <int D, class T> class System;
      template <int D, class T> class Simulator;
      template <int D, class T> class BdSimulator;
   }
}


namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Template for abstract base class for Brownian dynamics (BD) steps.
   *
   * The virtual step() method must generate a single BD step. The step
   * generally includes a random change in field components other than
   * the pressure-like field, followed by use of a Compressor to adjust
   * the pressure-like field to re-establish a homogeneous total monomer
   * concentration.
   *
   * Template parameters:
   *
   *   - D : dimension of space (D=1, 2, or 3)
   *   - T : backend identifier class (CPT or CUT)
   *
   * \see \ref psfts_algo_brownian_page "Manual Page"
   * \ingroup Rp_Fts_Brownian_Module
   */
   template <int D, class T>
   class BdStep : public ParamComposite
   {
   public:

      /**
      * Constructor.
      *
      * \param simulator  parent BdSimulator object
      */
      BdStep(BdSimulator<D,T>& simulator);

      ~BdStep() = default;

      // Prohibit copying and assignment
      BdStep(BdStep<D,T> const &) = delete;
      BdStep<D,T>& operator = (BdStep<D,T> const &) = delete;

      /**
      * Read required parameters from file.
      *
      * Empty default implementation.
      *
      * \param in  input parameter stream
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
      * The default implementation returns true.
      *
      * \return true to save, or false otherwise
      */
      virtual bool needsDc()
      {  return true; }

      /**
      * Output any timing results to ostream.
      *
      * Default implementation does nothing.
      *
      * \param out output stream
      */
      virtual void outputTimers(std::ostream& out);

      /**
      * Clear timers.
      *
      * Default implementation does nothing.
      */
      virtual void clearTimers();

      /**
      * Output statistics for this stepper (after end of simulation).
      *
      * Default implementation does nothing.
      */
      virtual void output();

   protected:

      /**
      * Get parent System object.
      */
      System<D,T>& system();

      /**
      * Get parent BdSimulator object.
      */
      BdSimulator<D,T>& simulator();

      /**
      * Get scalar random number generator of parent simulator.
      */
      Random& random();

      /**
      * Get vector random number generator of parent simulator.
      */
      typename T::VecRandom& vecRandom();

   private:

      /// Pointer to parent BdSimulator object
      BdSimulator<D,T>* simulatorPtr_;

      /// Pointer to parent System object
      System<D,T>* systemPtr_;

      /// Pointer to the scalar random number generator
      Random  *randomPtr_;

      /// Pointer to the vector random number generator
      typename T::VecRandom  *vecRandomPtr_;

   };

   // Protected inline methods

   /*
   * Get parent System object.
   */
   template <int D, class T> inline
   System<D,T>& BdStep<D,T>::system()
   {  return *systemPtr_; }

   /*
   * Get parent BdSimulator object.
   */
   template <int D, class T> inline
   BdSimulator<D,T>& BdStep<D,T>::simulator()
   {  return *simulatorPtr_; }

   /*
   * Get the scalar random number generator.
   */
   template <int D, class T> inline
   Random& BdStep<D,T>::random()
   {  return *randomPtr_; }

   /*
   * Get the vector random number generator.
   */
   template <int D, class T> inline
   typename T::VecRandom& BdStep<D,T>::vecRandom()
   {  return *vecRandomPtr_; }

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(BdStep)
}
}
#endif
