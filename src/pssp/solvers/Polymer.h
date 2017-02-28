#ifndef PSSP_POLYMER_H
#define PSSP_POLYMER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Block.h"
#include <pscf/solvers/PolymerTmpl.h>

namespace Pscf { 
namespace Pssp { 

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
   * \ingroup Pscf_Fd1d_Module
   */
   template <int D>
   class Polymer : public PolymerTmpl< Block<D> >
   {

   public:

      Polymer();

      ~Polymer();

      void setPhi(double phi);

      void setMu(double mu);

   };

}
}
#include "Polymer.tpp"
#endif
