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
   * The block concentrations stored in the constituent Block objects
   * contain the block concentrations (i.e., volume fraction) fields
   * computed in the most recent call of the compute function.
   *
   * The phi() and mu() accessor functions, which are inherited from
   * Pscf::PolymerSpecies, return the value of phi (spatial average volume
   * fraction of a species) or mu (species chemical potential) computed in 
   * the last call of the compute function.  If the ensemble for this 
   * species is closed, phi is read from the parameter file and mu is 
   * computed. If the ensemble is open, mu is read from the parameter 
   * file and phi is computed.
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
      * are computed for all propagators and blocks, along with molecular
      * partition function q and phi or mu.
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
