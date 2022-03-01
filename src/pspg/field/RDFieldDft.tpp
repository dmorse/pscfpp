#ifndef PSPG_R_DFIELD_DFT_TPP
#define PSPG_R_DFIELD_DFT_TPP

/*
* PSCF++ Package 
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "RDFieldDft.h"

namespace Pscf {
namespace Pspg {

   using namespace Util;

   /**
   * Default constructor.
   */
   template <int D>
   RDFieldDft<D>::RDFieldDft()
    : DField<cudaComplex>()
   {}

   /*
   * Destructor.
   */
   template <int D>
   RDFieldDft<D>::~RDFieldDft()
   {}

   /*
   * Copy constructor.
   *
   * Allocates new memory and copies all elements by value.
   *
   *\param other the RField<D> to be copied.
   */
   template <int D>
   RDFieldDft<D>::RDFieldDft(const RDFieldDft<D>& other)
    : DField<cudaComplex>(other)
   {
      meshDimensions_ = other.meshDimensions_;
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
   RDFieldDft<D>& RDFieldDft<D>::operator = (const RDFieldDft<D>& other)
   {
      
      DField<cudaComplex>::operator = (other);
      meshDimensions_ = other.meshDimensions_;

      return *this;
   }

}
}
#endif
