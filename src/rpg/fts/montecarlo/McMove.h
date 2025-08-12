#ifndef RPG_MC_MOVE_H
#define RPG_MC_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include <util/param/ParamComposite.h>
#include <util/random/Random.h>
#include <pscf/cuda/CudaRandom.h>
#include <util/global.h>
#include <util/misc/Timer.h>

namespace Pscf {
namespace Rpg
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
   * \ingroup Rpg_Fts_MonteCarlo_Module
   */
   template <int D>
   class McMove : public ParamComposite
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent McSimulator object
      */
      McMove(McSimulator<D>& simulator);

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
      * \return true if converged, false if failed to converge.
      */
      virtual bool move();
      
      /**
      * Decide whether cc fields need to be saved for move
      * The default implementation is false
      */
      virtual bool needsCc()
      {  return false; }
      
      /**
      * Decide whether dc fields need to be saved for move
      * The default implementation is false
      */
      virtual bool needsDc()
      { return false; }
      
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
      * Return number of moves that fail to converge.
      */
      long nFail() const;

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
      * Increment the number of failed moves.
      */
      void incrementNFail();

      /**
      * Get parent System object.
      */
      System<D>& system();

      /**
      * Get parent McSimulator object.
      */
      McSimulator<D>& simulator();

      /**
      * Get Random number generator of parent System.
      */
      Random& random();
      
      /**
      * Get cuda random number generator by reference.
      */
      CudaRandom& cudaRandom();
      
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
      Timer componentTimer_;
      Timer hamiltonianTimer_;
      Timer decisionTimer_;
      Timer totalTimer_;

   private:

      /// Pointer to parent McSimulator object
      McSimulator<D>* simulatorPtr_;

      /// Pointer to parent System object
      System<D>* systemPtr_;

      /// Pointer to random number generator
      Random  *randomPtr_;
      
      /// Pointer to cudaRandom number generator
      CudaRandom  *cudaRandomPtr_;

      /// Probability of choosing this move
      double  probability_;

      /// Number of moves that have been attempted by this object.
      long  nAttempt_;

      /// Number of moves that have been accepted by this object.
      long  nAccept_;
      
      /// Number of moves that fail to converge.
      long  nFail_;
      
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
   
   /*
   * Return number of moves that fail to converge.
   */
   template <int D>
   inline long McMove<D>::nFail() const
   {  return nFail_; }

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
   * Increment the number of fail moves.
   */
   template <int D>
   inline void McMove<D>::incrementNFail()
   {  ++nFail_; }

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
   inline McSimulator<D>& McMove<D>::simulator()
   {  return *simulatorPtr_; }

   /*
   * Get Random number generator.
   */
   template <int D>
   inline Random& McMove<D>::random()
   {  return *randomPtr_; }

   /*
   * Get CudaRandom number generator.
   */
   template <int D>
   inline CudaRandom& McMove<D>::cudaRandom()
   {  return *cudaRandomPtr_; }
   
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

   #ifndef RPG_MC_MOVE_TPP
   // Suppress implicit instantiation
   extern template class McMove<1>;
   extern template class McMove<2>;
   extern template class McMove<3>;
   #endif

}
}
#endif
