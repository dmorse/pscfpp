#ifndef PRDC_CPU_FIELD_BASIS_CONVERTER_TPP
#define PRDC_CPU_FIELD_BASIS_CONVERTER_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldBasisConverter.h"
#include <prdc/cpu/RField.h>
#include <prdc/cpu/RFieldDft.h>
#include <util/containers/DArray.h>

namespace Pscf {
namespace Prdc {
namespace Cpu {

   using namespace Util;

   /*
   * Constructor.
   */
   template <int D>
   FieldBasisConverter<D>::FieldBasisConverter(DMatrix<double> basis,
                                               double normSq)
    : basis_(basis),
      normSq_(normSq),
      nMonomer_(basis.capacity1())
   {
      UTIL_CHECK(basis.isAllocated());
      UTIL_CHECK(basis_.capacity1() == basis_.capacity2());
      UTIL_CHECK(normSq_ > 0.0);
      UTIL_CHECK(nMonomer_ > 1);
   }

   /*
   * Constructor.
   */
   template <int D>
   FieldBasisConverter<D>::FieldBasisConverter(DMatrix<double> basis,
                                               int normSq)
    : basis_(basis),
      normSq_(double(normSq)),
      nMonomer_(basis.capacity1())
   {
      UTIL_CHECK(basis.isAllocated());
      UTIL_CHECK(basis_.capacity1() == basis_.capacity2());
      UTIL_CHECK(normSq_ > 0.0);
      UTIL_CHECK(nMonomer_ > 1);
   }

   /*
   * Destructor.
   */
   template <int D>
   FieldBasisConverter<D>::~FieldBasisConverter()
   {}

   /*
   * Compute pointwise components of r-grid fields in basis.
   */
   template <int D> 
   void 
   FieldBasisConverter<D>::convertToBasis(DArray<RField<D>> const & in,
                                          DArray<RField<D>> & out) const
   {
      // Preconditions
      UTIL_CHECK(in.isAllocated());
      UTIL_CHECK(out.isAllocated());
      UTIL_CHECK(in.capacity() == nMonomer_);
      UTIL_CHECK(out.capacity() == nMonomer_);
      int meshSize = in[0].capacity();
      for (int i=0; i < nMonomer_; ++i) {
         UTIL_CHECK(in[i].capacity() == meshSize);
         UTIL_CHECK(out[i].capacity() == meshSize);
      }

      int i, j, k;

      // Loop over components in basis
      for (i = 0; i < nMonomer_; ++i) {
         RField<D> wo = out[i];

         // Initialize wo = out[i] to zero
         for (k = 0; k < meshSize; ++k) {
            wo[k] = 0.0; 
         }

         // Loop over monomer types
         for (j = 0; j < nMonomer_; ++j) {
            double vec = basis_(i, j)/normSq_;
            RField<D> wi = in[j];

            // Loop over grid points
            for (k = 0; k < meshSize; ++k) {
               wo[k] += vec*wi[k];
            }

         } // for j < nMonomer_
      } // for i < nMonomer_

   }

   /*
   * Compute r-grid fields from components in basis.
   */
   template <int D> 
   void 
   FieldBasisConverter<D>::convertFromBasis(DArray<RField<D>> const & in,
                                            DArray<RField<D>> & out) const
   {
      // Preconditions
      UTIL_CHECK(in.isAllocated());
      UTIL_CHECK(out.isAllocated());
      UTIL_CHECK(in.capacity() == nMonomer_);
      UTIL_CHECK(out.capacity() == nMonomer_);
      int meshSize = in[0].capacity();
      for (int i=0; i < nMonomer_; ++i) {
         UTIL_CHECK(in[i].capacity() == meshSize);
         UTIL_CHECK(out[i].capacity() == meshSize);
      }

      int i, j, k;
      // Loop over monomer types
      for (i = 0; i < nMonomer_; ++i) {
         RField<D> wo = out[i];

         // Initialize wo = out[i] to zero
         for (k = 0; k < meshSize; ++k) {
            wo[k] = 0.0; 
         }

         // Loop over components in basis
         for (j = 0; j < nMonomer_; ++j) {
            double vec = basis_(j, i);
            RField<D> wi = in[j];

            // Loop over grid points
            for (k = 0; k < meshSize; ++k) {
               wo[k] += vec*wi[k];
            }

         } // for j < nMonomer_
      } // for i < nMonomer_

   }

}
}
}
#endif
