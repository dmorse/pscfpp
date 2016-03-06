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
      UTIL_CHECK(nx_ == f.capacity());

      GeometryMode mode = geometryMode();
      double sum = 0.0;
      double norm = 0.0;
      if (mode == Planar) {

         sum += 0.5*f[0];
         for (int i = 1; i < nx_ - 1; ++i) {
            sum += f[i];
         }
         sum += 0.5*f[nx_ - 1];
         norm = double(nx_ - 1);

      } else 
      if (mode == Cylindrical) {

         // First value
         double x0 = xMin_/dx_;
         if (x0 < 0.1) {
            sum  += f[0]/8.0;
            norm += 1.0/8.0;
         } else {
            sum += 0.5*x0*f[0];
            norm += 0.5*x0;
         }

         // Interior values
         double x;
         for (int i = 1; i < nx_ - 1; ++i) {
            x = x0 + double(i);
            sum  += x*f[i];
            norm += x;
         }

         // Last value
         x = x0 + double(nx_-1);
         sum += 0.5*x*f[nx_-1];
         norm += 0.5*x;

      } else
      if (mode == Spherical) {

         // First value
         double x0 = xMin_/dx_;
         if (x0 < 0.1) {
            sum  += f[0]/24.0;
            norm += 1.0/24.0;
         } else {
            sum += 0.5*x0*x0*f[0];
            norm += 0.5*x0*x0;
         }

         // Interior values
         double x;
         for (int i = 1; i < nx_ - 1; ++i) {
            x = x0 + double(i);
            sum  += x*x*f[i];
            norm += x*x;
         }

         // Last value
         x = x0 + double(nx_-1);
         sum += 0.5*x*x*f[nx_-1];
         norm += 0.5*x*x;

      } else {

         UTIL_THROW("Invalid geometry mode");

      }

      return sum/norm;
   }
 
   /*
   * Compute inner product of two real fields.
   */
   double Domain::innerProduct(Field const & f, Field const & g) const
   {
      // Preconditions
      UTIL_CHECK(nx_ > 1);
      UTIL_CHECK(nx_ == f.capacity());
      UTIL_CHECK(nx_ == g.capacity());
      if (!work_.isAllocated()) {
         work_.allocate(nx_);
      } else {
         UTIL_ASSERT(nx_ == work_.capacity());
      } 

      // Compute average of f(x)*g(x)
      for (int i = 0; i < nx_; ++i) {
         work_[i] = f[i]*g[i];
      }
      return spatialAverage(work_);
   }

}
}
