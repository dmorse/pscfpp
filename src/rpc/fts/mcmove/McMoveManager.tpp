#ifndef RPC_MC_MOVE_MANAGER_TPP
#define RPC_MC_MOVE_MANAGER_TPP

/*
* PSCF - Polymer Self-Consistent Field 
*
* Copyright 2016 - 2023, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <rpc/fts/mcmove/McMoveManager.h>
#include <rpc/fts/mcmove/McMoveFactory.h>
#include <rpc/fts/mcmove/McSimulator.h>

#include <util/random/Random.h>

namespace Pscf {
namespace Rpc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   McMoveManager<D>::McMoveManager(McSimulator<D>& simulator, 
                                   System<D>& system)
   : Manager< McMove<D> >(),
     simulatorPtr_(&simulator),
     systemPtr_(&system),
     randomPtr_(&simulator.random())
   {  setClassName("McMoveManager"); }

   // Destructor
   template <int D>
   McMoveManager<D>::~McMoveManager()
   {}

   /*
   * Return a pointer to a new McMoveFactory object.
   */
   template <int D>
   Factory< McMove<D> >* McMoveManager<D>::newDefaultFactory() const
   {  return new McMoveFactory<D>(*simulatorPtr_); }

   /* 
   * Read instructions for creating objects from file.
   */
   template <int D>
   void McMoveManager<D>::readParameters(std::istream &in)
   {
      // Read parameters for all McMove<D> objects
      Manager< McMove<D> >::readParameters(in);

      // Allocate and store probabilities
      probabilities_.allocate(size());
      double  totalProbability = 0.0;
      int     iMove;
      for (iMove = 0; iMove < size(); ++iMove) {
         probabilities_[iMove] = (*this)[iMove].probability();
         totalProbability += probabilities_[iMove];
      }

      // Allocate and store and normalize probabilities
      for (iMove = 0; iMove < size(); ++iMove) {
         probabilities_[iMove] = probabilities_[iMove]/totalProbability;
         (*this)[iMove].setProbability(probabilities_[iMove]);
      }
   }

   /*
   * Initialize all moves just prior to a run.
   */
   template <int D>
   void McMoveManager<D>::setup()
   {
      for (int iMove = 0; iMove < size(); ++iMove) {
         (*this)[iMove].setup();
      }
   }

   /*
   * Choose a McMove at random.
   */
   template <int D>
   McMove<D>& McMoveManager<D>::chooseMove()
   {
      int iMove;
      iMove = randomPtr_->drawFrom(&probabilities_[0], size());
      return (*this)[iMove];
   }

   /*
   * Output statistics for every move.
   */
   template <int D>
   void McMoveManager<D>::output() const
   {
      for (int i=0; i< size(); i++) {
         (*this)[i].output();
      }
   }
   
   /*
   * Log output timing results 
   */
   template <int D>
   void McMoveManager<D>::outputTimers(std::ostream& out) const
   {
      for (int i=0; i< size(); i++) {
         (*this)[i].outputTimers(out);
      }
   }
   
   /*
   * Clear timers 
   */
   template <int D>
   void McMoveManager<D>::clearTimers()
   {
      for (int i=0; i< size(); i++) {
         (*this)[i].clearTimers();
      }
   }
   
   /*
   * Decide whether any move needs to store cc fields.
   */
   template <int D>
   bool McMoveManager<D>::needsCc()
   {
      for (int i=0; i< size(); i++) {
         if((*this)[i].needsCc()){
            return true;
         }
      }
      return false;
   }
   
   /*
   * Decide whether any move needs to store dc fields.
   */
   template <int D>
   bool McMoveManager<D>::needsDc()
   {
      for (int i=0; i< size(); i++) {
         if((*this)[i].needsDc()){
            return true;
         }
      }
      return false;
   }
   
}
}
#endif
