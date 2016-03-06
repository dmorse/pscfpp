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

   /*
   * Compute spatial average of a field.
   */
   double Domain::spatialAverage(Field const & f) const
   {
      UTIL_CHECK(nx_ > 1);
      UTIL_CHECK(dx_ > 0.0);
      UTIL_CHECK(xMax_ - xMin_ >=  dx_);
      UTIL_CHECK(f.capacity() == nx_);

      double sum = 0.0;
      double norm = 0.0;
      if (geometryMode() == Planar) {
         for (int i = 1; i < nx_ - 1; ++i) {
            sum += f[i];
         }
         #ifdef FD1D_DOMAIN_IDENTITY_NORM
         sum += f[0];
         sum += f[nx_ - 1];
         norm = double(nx_);
         #else
         sum += 0.5*f[0];
         sum += 0.5*f[nx_ - 1];
         norm = double(nx_ - 1);
         #endif
      } else 
      if (geometryMode() == Cylindrical) {
         double x0 = xMin_/dx_;
         sum += 0.5*x0*f[0];
         double x;
         for (int i = 1; i < nx_ - 1; ++i) {
            x = x0 + double(i);
            sum  += x*f[i];
            norm += x;
         }
         x = x0 + double(nx_-1);
         sum += 0.5*x*f[nx_-1];
         norm += 0.5*x;
      } else
      if (geometryMode() == Spherical) {
         UTIL_THROW("Spherical average not yet implemented");
      }
      return sum/norm;
   }
 
   /*
   * Compute inner product of two real fields.
   */
   double Domain::innerProduct(Field const & f, Field const & g) const
   {
      double sum = 0.0;
      double norm = 0.0;
      if (geometryMode() == Planar) {
         for (int i = 1; i < nx_ - 1; ++i) {
            sum += f[i]*g[i];
         }
         #ifdef FD1D_DOMAIN_IDENTITY_NORM
         sum += f[0]*g[0];
         sum += f[nx_ - 1]*g[nx_ - 1];
         norm = double(nx_);
         #else
         sum += 0.5*f[0]*g[0];
         sum += 0.5*f[nx_ - 1]*g[nx_ - 1];
         norm = double(nx_ - 1);
         #endif
      } else {
         UTIL_THROW("Non-planar inner product not yet implemented");
      }
      return sum/norm;
   }

}
}
