#ifndef PSPC_MC_MOVE_MANAGER_TPP
#define PSPC_MC_MOVE_MANAGER_TPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <pspc/mcmove/McMoveManager.h>
#include <pspc/mcmove/McMoveFactory.h>
#include <pspc/mcmove/McSimulator.h>

#include <util/random/Random.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   McMoveManager<D>::McMoveManager(McSimulator<D>& mcSimulator, 
                                   System<D>& system)
   : Manager< McMove<D> >(),
     mcSimulatorPtr_(&mcSimulator),
     systemPtr_(&system),
     randomPtr_(&mcSimulator.random())
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
   {  return new McMoveFactory<D>(*mcSimulatorPtr_); }

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
   void McMoveManager<D>::output()
   {
      for (int i=0; i< size(); i++) {
         (*this)[i].output();
      }
   }

}
}
#endif
