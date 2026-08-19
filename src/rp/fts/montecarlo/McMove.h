#ifndef RP_MC_MOVE_H
#define RP_MC_MOVE_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>  // base class
#include <util/misc/Timer.h>            // member
#include <util/global.h>
#include <pscf/backends/TmplDeclare.h>

// Forward declarations
namespace Util { 
   class Random; 
}
namespace Pscf {
   namespace Rp {
      template <int D, class T> class System;
      template <int D, class T> class Simulator;
      template <int D, class T> class McSimulator;
   }
}


namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * McMove is an abstract base class for Monte Carlo moves.
   *
   * The virtual move() function must generate a trial move, decide whether
   * to accept or reject it, and update the associated System fields if
   * the move is accepted.
   *
   * Specializations of this class template are used as base classes for 
   * two closely analogous class templates, both also named McMove, that 
   * are defined in Rpc and Rpg namespaces and used in the pscf_rpc and 
   * pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *    - D : dimension
   *    - T : backend identifier class (CPT or CUT)
   *
   * \ingroup Rp_Fts_MonteCarlo_Module
   */
   template <int D, class T>
   class McMove : public ParamComposite
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator  parent McSimulator object
      */
      McMove(McSimulator<D,T>& simulator);

      /**
      * Destructor.
      */
      ~McMove() = default;

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
      * This method zeros the statistical accumulators nAttempt, nAccept,
      * and nFail. Derived class re-implementations should call this
      * default version and then complete any other initialization that
      * is required at the beginning of each simulation run within a
      * possible sequence of several such runs.
      */
      virtual void setup();

      /**
      * Generate, attempt, and accept or reject a Monte Carlo move.
      *
      * Implementations of this by subclasses should:
      *
      *    - Generate a proposed move
      *    - Decide to accept or reject (using Random::metropolis)
      *    - Restore the old system state if rejected, or update
      *      the system state if accepted.
      *
      * The default implementation provides a skeleton that calls
      * the virtual attemptMove() function, calls the compressor
      * after attemptMove(), and uses a Metropolis test based on
      * the change of Hamiltonian for acceptance or rejection. 
      * MC moves for which this skeleton is appropriate can be
      * implemented by redefining the attemptMove() function and
      * using it in this default implementation of move(). MC moves
      * for which this skeleton is inappropriate or inadequate can
      * be implemented by redefining the move function.
      *
      * \return true if converged, false if failed to converge.
      */
      virtual bool move();

      /**
      * Decide whether cc fields need to be saved for move.
      *
      * The default implementation returns false.
      */
      virtual bool needsCc()
      {  return false; }

      /**
      * Decide whether dc fields need to be saved for move.
      *
      * The default implementation returns false.
      */
      virtual bool needsDc()
      {  return false; }

      /**
      * Write timing results to a file.
      *
      * \param out  output stream
      */
      virtual void outputTimers(std::ostream& out);

      /**
      * Clear timers
      */
      virtual void clearTimers();

      // Accessor Functions

      /**
      * Return the probability of choosing this McMove.
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
      * Return number of moves that failed to converge.
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
      * Write timing data to a file, without a class label.
      *
      * \param out  output stream
      */
      void outputTimerData(std::ostream& out);

      /**
      * Get parent System object (non-const ref).
      */
      System<D,T>& system();

      /**
      * Get parent System object (const ref)
      */
      System<D,T> const & system() const;

      /**
      * Get parent McSimulator object (non-const ref).
      */
      McSimulator<D,T>& simulator();

      /**
      * Get parent McSimulator object (const ref)
      */
      McSimulator<D,T> const & simulator() const;

      /**
      * Get the scalar random number generator.
      */
      Random& random();

      /**
      * Get the vector random number generator.
      */
      typename T::VecRandom& vecRandom();

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

      /// Pointer to parent McSimulator object.
      McSimulator<D,T>* simulatorPtr_;

      /// Pointer to parent System object.
      System<D,T>* systemPtr_;

      /// Pointer to the scalar random number generator.
      Random* randomPtr_;

      /// Pointer to the vector random number generator.
      typename T::VecRandom* vecRandomPtr_;

      /// Probability of choosing this move
      double probability_;

      /// Number of moves that have been attempted by this object.
      long nAttempt_;

      /// Number of moves that have been accepted by this object.
      long nAccept_;

      /// Number of moves that fail to converge.
      long nFail_;

   };

   // Public inline methods

   /*
   * Return number of moves that have been attempted.
   */
   template <int D, class T> inline
   long McMove<D,T>::nAttempt() const
   {  return nAttempt_; }

   /*
   * Return number of moves that have been accepted.
   */
   template <int D, class T> inline
   long McMove<D,T>::nAccept() const
   {  return nAccept_; }

   /*
   * Return number of moves that fail to converge.
   */
   template <int D, class T> inline
   long McMove<D,T>::nFail() const
   {  return nFail_; }

   // Protected inline methods

   /*
   * Increment the number of attempted moves.
   */
   template <int D, class T>
   inline void McMove<D,T>::incrementNAttempt()
   {  ++nAttempt_; }

   /*
   * Increment the number of accepted moves.
   */
   template <int D, class T>
   inline void McMove<D,T>::incrementNAccept()
   {  ++nAccept_; }

   /*
   * Increment the number of fail moves.
   */
   template <int D, class T>
   inline void McMove<D,T>::incrementNFail()
   {  ++nFail_; }

   /*
   * Get parent System object (non-const ref).
   */
   template <int D, class T> inline
   System<D,T>& McMove<D,T>::system()
   {  return *systemPtr_; }

   /*
   * Get parent System object (const ref).
   */
   template <int D, class T> inline
   System<D,T> const & McMove<D,T>::system() const
   {  return *systemPtr_; }

   /*
   * Get parent McSimulator object (non-const ref).
   */
   template <int D, class T> inline
   McSimulator<D,T>& McMove<D,T>::simulator()
   {  return *simulatorPtr_; }

   /*
   * Get parent McSimulator object (const ref).
   */
   template <int D, class T> inline
   McSimulator<D,T> const & McMove<D,T>::simulator() const
   {  return *simulatorPtr_; }

   /*
   * Get the scalar andom number generator.
   */
   template <int D, class T> inline
   Random& McMove<D,T>::random()
   {  return *randomPtr_; }

   /*
   * Get the vector random number generator.
   */
   template <int D, class T> inline
   typename T::VecRandom& McMove<D,T>::vecRandom()
   {  return *vecRandomPtr_; }

   /*
   * Get the probability.
   */
   template <int D, class T> inline
   double McMove<D,T>::probability() const
   {  return probability_; }

   /*
   * Set the probability.
   */
   template <int D, class T> inline
   void McMove<D,T>::setProbability(double probability)
   {  probability_ = probability; }

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE(McMove)

}
}
#endif
