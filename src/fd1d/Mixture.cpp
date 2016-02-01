/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"
#include "Grid.h"

#include <cmath>

namespace Pscf { 
namespace Fd1d
{ 

   Mixture::Mixture()
   {  setClassName("Mixture"); }

   Mixture::~Mixture()
   {}

   void Mixture::readParameters(std::istream& in)
   {
      MixtureTmpl<Polymer, Solvent>::readParameters(in);
      read(in, "ds", ds_);
   }

   void Mixture::setGrid(Grid const& grid)
   {
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer()+ nSolvent() > 0);
      UTIL_CHECK(ds_ > 0);

      gridPtr_ = &grid;
      double nx = grid.nx();
      int i, j;

      // Allocate chemical potential and concentration fields
      for (int i = 0; i < nMonomer(); ++i) {
         wField(i).allocate(nx);
         cField(i).allocate(nx);
      }

      // Set discretization for all blocks
      for (i = 0; i < nPolymer(); ++i) {
         for (j = 0; j < polymer(i).nBlock(); ++j) {
            polymer(i).block(j).setDiscretization(grid, ds_);
         }
      }

   }

   void Mixture::compute()
   {
      UTIL_CHECK(gridPtr_);
      UTIL_CHECK(grid().nx() > 0);
      UTIL_CHECK(nMonomer() > 0);
      UTIL_CHECK(nPolymer() + nSolvent() > 0);

      int nx = grid().nx();
      int i, j, k;

      // Clear all monomer concentration fields
      for (i = 0; i < nMonomer(); ++i) {
         for (j = 0; j < nx; ++j) {
            cField(i)[j] = 0.0;
         }
      }

      // Solve MDE for all polymers
      for (i = 0; i < nPolymer(); ++i) {
         polymer(i).compute(wFields());
      }

      // Accumulate monomer concentration fields
      for (i = 0; i < nPolymer(); ++i) {
         for (j = 0; j < polymer(i).nBlock(); ++j) {
            int monomerId = polymer(i).block(j).monomerId();
            UTIL_CHECK(monomerId >= 0);
            UTIL_CHECK(monomerId < nMonomer());
            CField& monomerField = cField(monomerId);
            CField& blockField = polymer(i).block(j).cField();
            for (k = 0; k < nx; ++k) {
               monomerField[k] += blockField[k];
            }
         }
      }
    
   }
}
}
