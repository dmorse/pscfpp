#ifndef RPC_MC_MOVE_H
#define RPC_MC_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include <util/param/ParamComposite.h>
#include <util/random/Random.h>
#include <util/global.h>
#include <util/misc/Timer.h>

namespace Pscf {
namespace Rpc
{

   using namespace Util;

   template <int D> class System;
   template <int D> class McSimulator;

   /**
   * McMove is an abstract base class for Monte Carlo moves.
   *
   * The virtual move() method must generate a trial move, decide whether
   * to accept or reject it, and update the associated System fields
   * if it is accepted.
   *
   * \ingroup Rpc_Simulate_McMove_Module
   */
   template <int D>
   class McMove : public ParamComposite
   {

   public:

      /**
      * Constructor.
      *
      * \param mcSimulator  parent McSimulator object
      */
      McMove(McSimulator<D>& mcSimulator);

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      virtual ~McMove();

      /**
      * Read required parameters from file.
      *
      * Empty default implementation.
      */
      virtual void readParameters(std::istream &in);

      /**
      * Set the probability for this McMove.
      *
      * \param probability Probability of this move being chosen. 
      */
      void setProbability(double probability);

      /**
      * Setup before the beginning of each simulation run.
      *
      * This method zeros the statistical accumulators nAttempt
      * and nAccept. Derived class implementations should 
      * complete any other initialization that is required at 
      * the beginning of each simulation run within a sequence
      * of several such runs.
      */
      virtual void setup();

      /**
      * Generate, attempt, and accept or reject a Monte Carlo move.
      *
      * Implementations of this by subclasses should:
      *
      *    - Generate a proposed move 
      *    - Decide to accept or reject (use random::metropolis)
      *    - Restore the old system state if rejected, and update
      *      the system state if accepted. 
      *
      * The default implementation provides a skeleton that calls 
      * the virtual attemptMove() function, calls the compressor
      * after attemptMove(), and uses a Metropolis test based on
      * the change of Hamiltonian for acceptance or rejection. MC 
      * moves for which this skeleton is appropriate can be 
      * implemented by redefining the attemptMove() function and
      * using it in this default implementation of move(). MC moves 
      * for which this skeleton is inappropriate or inadequate can 
      * be implemented by redefining the move function. 
      *
      * \return true if accepted, false if rejected
      */
      virtual bool move();
      
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
      * Return probability for this McMove.
      */
      double probability() const;

      /**
      * Return number of moves that have been attempted.
      */
      long nAttempt() const;

      /**
      * Return number of moves that have been accepted.
      */
      long nAccept() const;

      /**
      * Output statistics for this move (at the end of simulation)
      */
      virtual void output();

   protected:

      /**
      * Increment the number of attempted moves.
      */
      void incrementNAttempt();

      /**
      * Increment the number of accepted moves.
      */
      void incrementNAccept();

      /**
      * Get parent System object.
      */
      System<D>& system();

      /**
      * Get parent McSimulator object.
      */
      McSimulator<D>& mcSimulator();

      /**
      * Get Random number generator of parent System.
      */
      Random& random();

      /**
      * Read the probability from file.
      */
      void readProbability(std::istream& in);

      /**
      *  Attempt unconstrained move.
      *
      *  This function should modify the system w fields in r-grid 
      *  format, as returned by system().w().rgrid(), in order apply 
      *  an unconstrained attempted move. The compressor will then be
      *  applied in order to restore the density constraint.
      *
      *  The default implementation is empty.
      */
      virtual void attemptMove() 
      {};

      /// Timers for McMove 
      Timer attemptMoveTimer_;
      Timer compressorTimer_;
      Timer computeWcTimer_;
      Timer computeHamiltonianTimer_;
      Timer decisionTimer_;
      Timer totalTimer_;

   private:

      /// Pointer to parent McSimulator object
      McSimulator<D>* mcSimulatorPtr_;

      /// Pointer to parent System object
      System<D>* systemPtr_;

      /// Pointer to random number generator
      Random  *randomPtr_;

      /// Probability of choosing this move
      double  probability_;

      /// Number of moves that have been attempted by this object.
      long  nAttempt_;

      /// Number of moves that have been accepted by this object.
      long  nAccept_;
      
   };

   // Public inline methods

   /*
   * Return number of moves that have been attempted.
   */
   template <int D>
   inline long McMove<D>::nAttempt() const
   {  return nAttempt_; }

   /*
   * Return number of moves that have been accepted.
   */
   template <int D>
   inline long McMove<D>::nAccept() const
   {  return nAccept_; }

   // Protected inline methods

   /*
   * Increment the number of attempted moves.
   */
   template <int D>
   inline void McMove<D>::incrementNAttempt()
   {  ++nAttempt_; }

   /*
   * Increment the number of accepted moves.
   */
   template <int D>
   inline void McMove<D>::incrementNAccept()
   {  ++nAccept_; }

   /*
   * Get parent System object.
   */
   template <int D>
   inline System<D>& McMove<D>::system()
   {  return *systemPtr_; }

   /*
   * Get parent McSimulator object.
   */
   template <int D>
   inline McSimulator<D>& McMove<D>::mcSimulator()
   {  return *mcSimulatorPtr_; }

   /*
   * Get Random number generator.
   */
   template <int D>
   inline Random& McMove<D>::random()
   {  return *randomPtr_; }

   /*
   * Get the probability.
   */
   template <int D>
   inline double McMove<D>::probability() const
   {  return probability_; }

   /*
   * Set the probability.
   */
   template <int D>
   inline void McMove<D>::setProbability(double probability)
   {  probability_ = probability; }

   #ifndef RPC_MC_MOVE_TPP
   // Suppress implicit instantiation
   extern template class McMove<1>;
   extern template class McMove<2>;
   extern template class McMove<3>;
   #endif

}
}
#endif
