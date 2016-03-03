/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Domain.h"

namespace Pscf { 
namespace Fd1d
{ 

   Domain::Domain()
    : xMin_(0.0),
      xMax_(0.0),
      dx_(0.0),
      nx_(0),
      geometryMode_(Planar)
   {  setClassName("Domain"); }

   Domain::~Domain()
   {}

   void Domain::setParameters(double xMin, double xMax, int nx)
   {
      xMin_ = xMin;
      xMax_ = xMax;
      nx_ = nx;
      dx_ = (xMax_ - xMin_)/double(nx_ - 1);
   }

   void Domain::readParameters(std::istream& in)
   {
      read(in, "xMin", xMin_);
      read(in, "xMax", xMax_);
      read(in, "nx", nx_);
      dx_ = (xMax_ - xMin_)/double(nx_ - 1);

      readOptional(in, "geometryMode", geometryMode_);
   }

}
}
