#ifndef PSPC_MC_MOVE_TPP
#define PSPC_MC_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "McMove.h"

#include <pspc/System.h>
#include <util/archives/Serializable_includes.h>

namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   McMove<D>::McMove(System<D>& system) 
    : systemPtr_(&system),
      randomPtr_(&system.random()),
      nAttempt_(0),
      nAccept_(0)
   {}

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   McMove<D>::~McMove()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void McMove<D>::readParameters(std::istream &in)
   {}
   
   /*
   * Read the probability from file.
   */
   template <int D>
   void McMove<D>::readProbability(std::istream &in)
   {  read<double>(in, "probability", probability_); }
  
   /*
   * Setup at beginning of loop.
   *
   * Trivial default implementation - initializes counters.
   */
   template <int D>
   void McMove<D>::setup()
   { 
      nAttempt_ = 0;
      nAccept_  = 0;
   }

   /*
   * Trivial default implementation - always returns false.
   */
   template <int D>
   bool McMove<D>::move()
   { 
      // Get free energy
      // Compute and store old Hamiltonian
      // Save old state
      // Attempt modification
      // Call compressor
      // Evaluate new Hamiltonian
      bool accept = false;
      // Call Random::metropolis
      // if (accept) {
      //    ++nAccept_;
      // } else {
      //    ++nReject
      //    Restore old state
      // }
      ++nAttempt_;
      return accept;
   }

   /*
   * Trivial default implementation - do nothing
   */
   template <int D>
   void McMove<D>::output()
   {}

}
}
#endif
