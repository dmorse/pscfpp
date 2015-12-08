#ifndef FD1D_PROPAGATOR_H
#define FD1D_PROPAGATOR_H

/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <chem/PropagatorTmpl.h>
#include <util/containers/DArray.h>

namespace Fd1d{ 

   using namespace Util;
   using namespace Chem;

   class Propagator : public PropagatorTmpl<Propagator>
   {

   public:

      // Public typedefs

      /**
      * Chemical potential field type.
      */
      typedef DArray<double> WField;

      /**
      * Monomer concentration field type.
      */
      typedef DArray<double> CField;

      // Member functions

      /**
      * Constructor.
      */
      Propagator();

      /**
      * Solve the modified diffusion equation for this block.
      *
      * \param w chemical potential field for appropriate monomer type
      */
      void solve(const WField& w);
  
      /**
      * Integrate to calculate monomer concentration for this block
      */ 
      void integrate(CField& integral);
   
   };

} 
#endif
