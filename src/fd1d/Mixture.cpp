/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Mixture.h"
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

      // Read discretization parameters
      read(in, "xMin", xMin_);
      read(in, "xMax", xMax_);
      read(in, "nx", nx_);
      read(in, "ds", ds_);
      dx_ = (xMax_ - xMin_)/double(nx_ - 1);

      // Allocate chemical potential and concentration fields
      for (int i = 0; i < nMonomer(); ++i) {
         wField(i).allocate(nx_);
         cField(i).allocate(nx_);
      } 

      // Set discretization for all blocks
      // double length;
      // int i, j, ns;
      int i, j;
      for (i = 0; i < nPolymer(); ++i) {
         for (j = 0; j < polymer(i).nBlock(); ++j) {
            #if 0
            length = polymer(i).block(j).length();
            UTIL_CHECK(length > 0);
            ns = floor(length/ds_ + 0.5) + 1;
            if (ns%2 == 0) {
               ns += 1;
            }
            #endif
            polymer(i).block(j).setDiscretization(xMin_, xMax_, nx_, ds_);
         }
      }

   }

   void Mixture::compute()
   {

      // Clear all monomer concentration fields
      int i, j;
      for (i = 0; i < nMonomer(); ++i) {
         for (j = 0; j < nx_; ++j) {
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
            for (int k = 0; k < nx_; ++k) {
               monomerField[k] += blockField[k];
            }
         }
      }
    
   }
}
}
