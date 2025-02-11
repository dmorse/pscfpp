/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Polymer.h"

namespace Pscf {
namespace R1d 
{ 

   Polymer::Polymer()
   {  setClassName("Polymer"); }

   Polymer::~Polymer()
   {}

   void Polymer::setPhi(double phi)
   {
      UTIL_CHECK(ensemble() == Species::Closed);  
      UTIL_CHECK(phi >= 0.0);  
      UTIL_CHECK(phi <= 1.0);  
      phi_ = phi; 
   }

   void Polymer::setMu(double mu)
   {
      UTIL_CHECK(ensemble() == Species::Open);  
      mu_ = mu; 
   }

   /*
   * Compute solution to MDE and concentrations.
   */ 
   void Polymer::compute(DArray<Block::WField> const & wFields)
   {

      // Setup solvers for all blocks
      int monomerId;
      for (int j = 0; j < nBlock(); ++j) {
         monomerId = block(j).monomerId();
         block(j).setupSolver(wFields[monomerId]);
      }

      // Solve MDE for all propagators
      solve();

      // Compute block concentration fields (thread model)
      double prefactor;
      prefactor = phi() / ( q() * length() );
      for (int i = 0; i < nBlock(); ++i) {
         block(i).computeConcentration(prefactor);
      }

   }

}
}
