/*
* PFTS - Polymer Field Theory Simulator
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "System.h"
#include <cmath>

namespace Fd1d
{ 

   System::System()
   {  setClassName("System"); }

   System::~System()
   {}

   void System::readParameters(std::istream& in)
   {
      SystemTmpl<Polymer, Solvent>::readParameters(in);

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

      // Set grid for all propagators
      double length;
      int i, j, ns;
      for (i = 0; i < nPolymer(); ++i) {
         for (j = 0; j < polymer(i).nBlock(); ++j) {
            length = polymer(i).block(j).length();
            UTIL_CHECK(length > 0);
            ns = floor(length/ds_ + 0.5) + 1;
            if (ns%2 == 0) {
               ns += 1;
            }
            polymer(i).propagator(j, 0).setGrid(xMin_, xMax_, nx_, ns);
            polymer(i).propagator(j, 1).setGrid(xMin_, xMax_, nx_, ns);
         }
      }

      // Allocate per-block concentration fields
      for (i = 0; i < nPolymer(); ++i) {
         for (j = 0; j < polymer(i).nBlock(); ++j) {
            polymer(i).blockCField(j).allocate(nx_);
         }
      }

   }

}
