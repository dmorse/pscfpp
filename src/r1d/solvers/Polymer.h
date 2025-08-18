#ifndef R1D_POLYMER_H
#define R1D_POLYMER_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <pscf/solvers/PolymerTmpl.h>   // Base class template
#include "Block.h"                      // Base class template parameter

namespace Pscf {
namespace R1d {

   using namespace Util;

   /**
   * Descriptor and solver for a block polymer species.
   *
   * A Polymer contains an array of Block objects, each of which contains
   * two propagator objects. After Polymer::compute is called, each block
   * stores the block concentration (volume fraction) fields computed for 
   * that block, and each Propagator contains the solution of the modified
   * diffusion equation for that propagator. Functions that return values
   * for molecular volume fraction phi, molecular chemical potential mu, 
   * and molecular partition function q are indirectly inherited from 
   * Pscf::PolymerSpecies.
   * 
   * Class R1d::Block is a subclass of Pscf::PolymerTmpl<Block>, which is
   * a subclass of Pscf::PolymerSpecies.
   *
   * \ref user_param_polymer_sec "Parameter File Format"
   * \ingroup R1d_Solver_Module
   */
   class Polymer : public PolymerTmpl<Block>
   {

   public:

      /**
      * Default constructor.
      */
      Polymer();

      /**
      * Destructor.
      */
      ~Polymer();

      /**
      * Compute solution to modified diffusion equation and concentrations.
      *
      * Upon return, propagator solutions and block concentration fields
      * have been computed for all propagators and blocks, along with the
      * molecular partition function q and phi or mu.
      *
      * \param wFields  array of chemica potential fields.
      */
      void compute(DArray< DArray<double> > const & wFields);

   private:

      // Restrict access
      using PolymerTmpl<Block>::solve;

   };

}
}
#endif
