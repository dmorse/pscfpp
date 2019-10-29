#ifndef FD1D_POLYMER_H
#define FD1D_POLYMER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.h"
#include <pscf/solvers/PolymerTmpl.h>

namespace Pscf { 
namespace Fd1d
{ 

   /**
   * Descriptor and solver for a branched polymer species.
   *
   * The block concentrations stored in the constituent Block
   * objects contain the block concentrations (i.e., volume 
   * fractions) computed in the most recent call of the compute 
   * function.
   *
   * The phi() and mu() accessor functions, which are inherited 
   * from PolymerTmp<Block>, return the value of phi (spatial 
   * average volume fraction of a species) or mu (chemical
   * potential) computed in the last call of the compute function.
   * If the ensemble for this species is closed, phi is read from 
   * the parameter file and mu is computed. If the ensemble is
   * open, mu is read from the parameter file and phi is computed.
   *
   * \ingroup Fd1d_Solver_Module
   */
   class Polymer : public PolymerTmpl<Block>
   {

   public:

      Polymer();

      ~Polymer();

      void setPhi(double phi);

      void setMu(double mu);

      /**
      * Compute solution to modified diffusion equation.
      *
      * Upon return, q functions and block concentration fields
      * are computed for all propagators and blocks. 
      */ 
      void compute(const DArray<Block::WField>& wFields);

   };

} 
}
#endif
