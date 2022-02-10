/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2019, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BFieldComparison.h"

namespace Pscf {
namespace Pspg {

   BFieldComparison::BFieldComparison(int begin)
    : FieldComparison< DArray<cudaReal*> > (begin)
   {}; 

   double BFieldComparison::compare(DArray<cudaReal*> const & a,
                                    DArray<cudaReal*> const & b, int nStar)
   {   
      UTIL_CHECK(a.capacity() > 0); 
      UTIL_CHECK(a.capacity() == b.capacity());
      UTIL_CHECK(nStar > 0);

      int m = a.capacity();
      double diff;
      maxDiff_ = 0.0;
      rmsDiff_ = 0.0;
      int i, j;
      for (i = 0; i < m; ++i) {
          
         for (j = begin_; j < nStar; ++j) {
            diff = std::abs(a[i][j] - b[i][j]);
            if (diff > maxDiff_) {
               maxDiff_ = diff;
            }   
            rmsDiff_ += diff*diff;
         }   
      }   
      rmsDiff_ = rmsDiff_/double(m*nStar);
      rmsDiff_ = sqrt(rmsDiff_);
      return maxDiff_;
   }  

}
}
