#ifndef PSPG_R_FIELD_DFT_TPP
#define PSPG_R_FIELD_DFT_TPP

/*
* PSCF Package 
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RFieldDft.h"

namespace Pscf {
namespace Pspg {

   using namespace Util;

   /**
   * Default constructor.
   */
   template <int D>
   RFieldDft<D>::RFieldDft()
    : Field<cudaComplex>()
   {}

   /*
   * Destructor.
   */
   template <int D>
   RFieldDft<D>::~RFieldDft()
   {}

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   *
   *\param other the RField<D> to be copied.
   */
   template <int D>
   RFieldDft<D>::RFieldDft(const RFieldDft<D>& other)
    : Field<cudaComplex>(other)
   {
      meshDimensions_ = other.meshDimensions_;
      dftDimensions_ = other.dftDimensions_;
   }

   /*
   * Assignment, element-by-element.
   *
   * This operator will allocate memory if not allocated previously.
   *
   * \throw Exception if other Field is not allocated.
   * \throw Exception if both Fields are allocated with unequal capacities.
   *
   * \param other the rhs Field
   */
   template <int D>
   RFieldDft<D>& RFieldDft<D>::operator = (const RFieldDft<D>& other)
   {
      
      Field<cudaComplex>::operator = (other);
      meshDimensions_ = other.meshDimensions_;
      dftDimensions_ = other.dftDimensions_;

      return *this;
   }

}
}
#endif
