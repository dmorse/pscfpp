/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Domain.h"
#include <util/math/Constants.h>

namespace Pscf { 
namespace Cyln
{

   using namespace Util; 

   Domain::Domain()
    : radius_(0.0),
      length_(0.0),
      volume_(0.0),
      dr_(0.0),
      dz_(0.0),
      nr_(0),
      nz_(0),
      nGrid_(0),
      work_()
   {  setClassName("Domain"); }

   Domain::~Domain()
   {}

   void Domain::readParameters(std::istream& in)
   {
      read(in, "radius", radius_);
      read(in, "length", length_);
      read(in, "nr", nr_);
      read(in, "nz", nz_);
      nGrid_ = nr_*nz_;
      dr_ = radius_/double(nr_ - 1);
      dz_ = length_/double(nz_ - 1);
      volume_ = length_*radius_*radius_*Constants::Pi;
   }

   void Domain::setParameters(double radius, double length, 
                              int nr, int nz)
   {
      UTIL_CHECK(radius > 0.0);
      UTIL_CHECK(length > 0.0);
      UTIL_CHECK(nr > 1);
      UTIL_CHECK(nz > 1);
      radius_ = radius;
      length_ = length;
      volume_ = length_*radius_*radius_*Constants::Pi;
      nr_ = nr;
      nz_ = nz;
      nGrid_ = nr_*nz_;
      dr_ = radius_/double(nr_ - 1);
      dz_ = length_/double(nz_ - 1);
   }

   /*
   * Compute spatial average of a field.
   */
   double Domain::spatialAverage(Field<double> const & f) const
   {
      UTIL_CHECK(nr_ == f.nr());
      UTIL_CHECK(nz_ == f.nz());
      UTIL_CHECK(nr_ > 1);
      UTIL_CHECK(nz_ > 1);
      UTIL_CHECK(dr_ > 0.0);
      UTIL_CHECK(dz_ > 0.0);
      UTIL_CHECK(radius_ >  dr_);

      double sum = 0.0;
      double norm = 0.0;
      double r;
      int i, j;

      int k = 0;
      // Central axis (r=0)
      r = 1.0/8.0;
      for (j = 0; j < nz_; ++j) {
         sum += f[k]*r;
         norm += r;
         ++k;
      }
      // Intermediate shells
      for (i = 1; i < nr_ - 1; ++i) {
         r = double(i);
         for (j = 0; j < nz_; ++j) {
            sum += r*f[0];
            norm += r;
            ++k;
         }
      }
      // Outer shell
      r = 0.5*double(nr_-1);
      for (j=0; j < nz_; ++j) {
         sum += r*f[nr_-1];
         norm += r;
         ++k;
      }
      UTIL_CHECK(k == nr_*nz_);

      return sum/norm;
   }
 
   /*
   * Compute inner product of two real fields.
   */
   double Domain::innerProduct(Field<double> const & f, Field<double> const & g) const
   {
      // Preconditions
      UTIL_CHECK(nGrid_ > 1);
      UTIL_CHECK(nGrid_ == f.capacity());
      UTIL_CHECK(nGrid_ == g.capacity());
      if (!work_.isAllocated()) {
         work_.allocate(nr_, nz_);
      } else {
         UTIL_ASSERT(nGrid_ == work_.capacity());
      } 

      // Compute average of f(x)*g(x)
      for (int i = 0; i < nGrid_; ++i) {
         work_[i] = f[i]*g[i];
      }
      return spatialAverage(work_);
   }

}
}
