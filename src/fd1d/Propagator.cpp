/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Propagator.h"

namespace Fd1d{

   using namespace Util;
   using namespace Chem;

   /*
   * Constructor.
   */
   Propagator::Propagator()
   {}

   /*
   * Solve the modified diffusion equation for this block.
   */
   void Propagator::solve(const Propagator::WField& w)
   {  
      setIsSolved(true); 
   }

   /*
   * Integrate to calculate monomer concentration for this block
   */
   void Propagator::integrate(Propagator::CField& integral)
   {
   }

}
