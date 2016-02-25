/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2013, David Morse (morse012@.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Grid.h"

namespace Pscf { 
namespace Fd1d
{ 

   Grid::Grid()
   {  setClassName("Grid"); }

   Grid::~Grid()
   {}

   void Grid::setParameters(double xMin, double xMax, int nx)
   {
      xMin_ = xMin;
      xMax_ = xMax;
      nx_ = nx;
      dx_ = (xMax_ - xMin_)/double(nx_ - 1);
   }

   void Grid::readParameters(std::istream& in)
   {
      read(in, "xMin", xMin_);
      read(in, "xMax", xMax_);
      read(in, "nx", nx_);
      dx_ = (xMax_ - xMin_)/double(nx_ - 1);
   }

}
}
