#ifndef PSPC_TEST_MOVE_TPP
#define PSPC_TEST_MOVE_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TestMove.h"
#include "McMove.h" 
#include <util/param/ParamComposite.h>
#include <pspc/System.h>
#include <util/archives/Serializable_includes.h>
#include <util/random/Random.h>


namespace Pscf {
namespace Pspc {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   TestMove<D>::TestMove(McSimulator<D>& simulator) 
    : McMove<D>(simulator)
   { setClassName("TestMove"); }

   /*
   * Destructor, empty default implementation.
   */
   template <int D>
   TestMove<D>::~TestMove()
   {}

   /*
   * ReadParameters, empty default implementation.
   */
   template <int D>
   void TestMove<D>::readParameters(std::istream &in)
   {
      //Read the probability
      readProbability(in);
      // attampt move range [A, -A]
      read(in, "A", A_);
   }
   
   /*
   * Attempt unconstrained move
   */
   template <int D>
   void TestMove<D>::attemptMove()
   {
      const int nMonomer = system().mixture().nMonomer();
      const int meshSize = system().domain().mesh().size();
      DArray< RField<D> > wField;
      wField.allocate(nMonomer);
      
      
      //New field is the w0 + the newGuess for the Lagrange multiplier field
      for (int i = 0; i < nMonomer; i++){
         wField[i].allocate(meshSize);
         for (int k = 0; k < meshSize; k++){
            //Random number generator
            double r = random().uniform(-A_,A_);
            wField[i][k] = system().w().rgrid()[i][k] + r;
         }
      }
      system().setWRGrid(wField);

   }


   /*
   * Trivial default implementation - do nothing
   */
   template <int D>
   void TestMove<D>::output()
   {}

}
}
#endif
