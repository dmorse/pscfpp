#ifndef CYLN_FIELD_TPP
#define CYLN_FIELD_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Field.h"

#include <fftw3.h>

namespace Pscf {
namespace Cyln
{

   using namespace Util;

   /*
   * Default constructor.
   */
   template <typename Data>
   Field<Data>::Field()
    : data_(0),
      capacity_(0),
      nr_(0),
      nz_(0),
      slices_()
   {}

   /*
   * Destructor.
   */
   template <typename Data>
   Field<Data>::~Field()
   {
      if (isAllocated()) {
         fftw_free(data_);
         capacity_ = 0;
      }
   }

   /*
   * Allocate the underlying C array.
   *
   * Throw an Exception if the Field has already allocated.
   *
   * \param capacity number of elements to allocate.
   */
   template <typename Data>
   void Field<Data>::allocate(int nr, int nz)
   {
      if (isAllocated()) {
         UTIL_THROW("Attempt to re-allocate a Field");
      }
      if (nr <= 0) {
         UTIL_THROW("Attempt to allocate with nr <= 0");
      }
      if (nz <= 0) {
         UTIL_THROW("Attempt to allocate with nz <= 0");
      }
      nr_ = nr;
      nz_ = nz;
      capacity_ = nr*nz;
      data_ = (Data*) fftw_malloc(sizeof(Data)*capacity_);

      // Allocate and associaetd array of 1D slices
      slices_.allocate(nr_);
      for (int i = 0; i < nr_; ++i) {
         slices_[i].associate(data_ + nz_*i, nz_);
      }

   }

   /*
   * Deallocate the underlying C array.
   *
   * Throw an Exception if this Field is not allocated.
   */
   template <typename Data>
   void Field<Data>::deallocate()
   {
      if (!isAllocated()) {
         UTIL_THROW("Array is not allocated");
      }
      fftw_free(data_);
      deallocate(slices_);
      capacity_ = 0;
      nr_ = 0;
      nz_ = 0;
   }

}
}
#endif
