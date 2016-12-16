/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Domain.h"
#include <util/math/Constants.h>

namespace Pscf { 
namespace Fd1d
{ 

   Domain::Domain()
    : xMin_(0.0),
      xMax_(0.0),
      dx_(0.0),
      volume_(0.0),
      nx_(0),
      mode_(Planar),
      isShell_(false)
   {  setClassName("Domain"); }

   Domain::~Domain()
   {}

   void Domain::readParameters(std::istream& in)
   {
      mode_ = Planar;
      read(in, "mode", mode_);
      if (mode_ != Planar) {
         read(in, "isShell", isShell_); 
      } 
      if (mode_ == Planar || isShell_) { 
         read(in, "xMin", xMin_);
      } else {
         xMin_ = 0.0;
      }
      read(in, "xMax", xMax_);
      read(in, "nx", nx_);
      dx_ = (xMax_ - xMin_)/double(nx_ - 1);
      computeVolume();
   }

   void Domain::setPlanarParameters(double xMin, double xMax, int nx)
   {
      mode_ = Planar;
      isShell_ = false;
      xMin_ = xMin;
      xMax_ = xMax;
      nx_ = nx;
      dx_ = (xMax_ - xMin_)/double(nx_ - 1);
      computeVolume();
   }

   void Domain::setShellParameters(GeometryMode mode, 
                                   double xMin, double xMax, int nx)
   {
      mode_ = mode;
      isShell_ = true;
      xMin_ = xMin;
      xMax_ = xMax;
      nx_ = nx;
      dx_ = (xMax_ - xMin_)/double(nx_ - 1);
      computeVolume();
   }

   void Domain::setCylinderParameters(double xMax, int nx)
   {
      mode_ = Cylindrical;
      isShell_ = false;
      xMin_ = 0.0;
      xMax_ = xMax;
      nx_ = nx;
      dx_ = xMax_/double(nx_ - 1);
      computeVolume();
   }

   void Domain::setSphereParameters(double xMax, int nx)
   {
      mode_ = Spherical;
      isShell_ = false;
      xMin_ = 0.0;
      xMax_ = xMax;
      nx_ = nx;
      dx_ = xMax_/double(nx_ - 1);
      computeVolume();
   }

   void Domain::computeVolume()
   {
      if (mode_ == Planar) {
         volume_ = xMax_ - xMin_;
      } else 
      if (mode_ == Cylindrical) {
         volume_ = xMax_*xMax_;
         if (isShell_) {
            volume_ -= xMin_*xMin_;
         }
         volume_ *= Constants::Pi;
      } else 
      if (mode_ == Spherical) {
         volume_ = xMax_*xMax_*xMax_;
         if (isShell_) {
            volume_ -= xMin_*xMin_*xMin_; 
         }
         volume_ *= Constants::Pi*(4.0/3.0);
      } else {
         UTIL_THROW("Invalid geometry mode");
      }
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

      double sum = 0.0;
      double norm = 0.0;
      if (mode_ == Planar) {

         sum += 0.5*f[0];
         for (int i = 1; i < nx_ - 1; ++i) {
            sum += f[i];
         }
         sum += 0.5*f[nx_ - 1];
         norm = double(nx_ - 1);

      } else 
      if (mode_ == Cylindrical) {

         double x0 = xMin_/dx_;

         // First value
         if (isShell_) {
            sum += 0.5*x0*f[0];
            norm += 0.5*x0;
         } else {
            sum += f[0]/8.0;
            norm += 1.0/8.0;
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
      if (mode_ == Spherical) {

         double x0 = xMin_/dx_;

         // First value
         if (isShell_) {
            sum += 0.5*x0*x0*f[0];
            norm += 0.5*x0*x0;
         } else {
            sum += f[0]/24.0;
            norm += 1.0/24.0;
         }

         // Interior values
         double x;
         for (int i = 1; i < nx_ - 1; ++i) {
            x = x0 + double(i);
            sum += x*x*f[i];
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
