#ifndef RP_MC_MOVE_MANAGER_H
#define RP_MC_MOVE_MANAGER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Manager.h>          // base class template
#include <rp/fts/montecarlo/McMove.h>    // template argument
#include <util/containers/DArray.h>      // member

// Forward declarations
namespace Util { 
   class Random; 
}
namespace Pscf {
   namespace Rp {
      template <int D, class T> class System;
      template <int D, class T> class Simulator;
      template <int D, class T> class McSimulator;
      template <int D, class T> class McMove;
   }
}


namespace Pscf {
namespace Rp {

   using namespace Util;

   /**
   * Manager for a set of McMove objects.
   *
   * Specializations of this class template are used as base classes for 
   * two closely analogous class templates, both also named McMoveManager,
   * that are defined in Rpc and Rpg namespaces for use in the pscf_rpc
   * and pscf_rpg programs, respectively.
   *
   * Template parameters:
   *
   *    - D : dimension
   *    - T : Types class, Rpc::Types<D> or Rpg::Types<D>
   *
   * \ingroup Rp_Fts_MonteCarlo_Module
   */
   template <int D, class T>
   class McMoveManager : public Manager< McMove<D,T> >
   {

   public:

      /**
      * Constructor.
      *
      * \param simulator parent McSimulator
      * \param system parent System
      */
      McMoveManager(McSimulator<D,T>& simulator, 
                    System<D,T>& system);

      /**
      * Destructor.
      */
      ~McMoveManager() = default;

      /**
      * Read instructions for creating McMove objects.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream &in);

      /**
      * Initialize at beginning of system run.
      *
      * This method calls the initialize method for every McMove.
      */
      void setup();

      /**
      * Choose an McMove at random, using specified probabilities.
      *
      * \return chosen McMove
      */
      McMove<D,T>& chooseMove();

      /**
      * Output statistics for all moves.
      */
      void output() const;

      /**
      * Return probability of move i.
      *
      * \param i index for McMove
      * \return probability of McMove number i
      */
      double probability(int i) const;

      /**
      * Write timing results to an output stream
      *
      * \param out  output stream (i.e., file) 
      */
      void outputTimers(std::ostream& out) const;

      /**
      * Clear timers.
      */
      void clearTimers();

      /**
      * Decide whether any move needs to store cc fields.
      */
      bool needsCc();

      /**
      * Decide whether any move needs to store dc fields.
      */
      bool needsDc();

   private:

      using McMoveT = McMove<D,T>;
      using Base = Manager< McMove<D,T> >;

      // Private data members

      /**
      * Array of McMove probabilities.
      */
      DArray<double>  probabilities_;

      /**
      * Pointer to parent Simulator.
      */
      McSimulator<D,T>* simulatorPtr_;

      /**
      * Pointer to parent System.
      */
      System<D,T>* systemPtr_;

      /**
      * Pointer to random number generator.
      */
      Random* randomPtr_;

      // Private member functions

      /**
      * Return pointer to a new McMoveFactory.
      */
      virtual Factory< McMove<D,T> >* newDefaultFactory() const;

   };

   // Inline functions

   /*
   * Return probability of move number i.
   */
   template <int D, class T> inline 
   double McMoveManager<D,T>::probability(int i) const
   {
      assert(i >= 0);
      assert(i < Base::size());
      return probabilities_[i];
   }

}
}
#endif
