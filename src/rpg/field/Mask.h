#ifndef RPG_MASK_H
#define RPG_MASK_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.h"              // parent class template parameter
#include <prdc/cuda/RField.h>     // parent class template parameter
#include <prdc/field/MaskTmpl.h>  // parent class

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /**
   * A field to which the total density is constrained.
   *
   * Please refer to the documentation of the base class Prdc::MaskTmpl
   * for more complete API documentation for this class template.
   * The public interface of Rpg::Mask is identical to that of the
   * base class template Prdc::MaskTmpl. 
   *
   * \ingroup Rpg_Field_Module
   */
   template <int D>
   class Mask : public Prdc::MaskTmpl< D, FieldIo<D>, Prdc::Cuda::RField<D> >
   {

   public:

      /**
      * Constructor.
      */
      Mask();

      /**
      * Destructor.
      */
      ~Mask();

   protected:

      /**
      * Calculate the average value of the rgrid_ member.
      */
      double rGridAverage() const;

   };

   #ifndef RPG_MASK_TPP
   extern template class Mask<1>;
   extern template class Mask<2>;
   extern template class Mask<3>;
   #endif

} // namespace Rpg

#ifndef RPG_MASK_TPP
namespace Prdc {
   extern template class MaskTmpl< 1, Rpg::FieldIo<1>, Cuda::RField<1> >;
   extern template class MaskTmpl< 2, Rpg::FieldIo<2>, Cuda::RField<2> >;
   extern template class MaskTmpl< 3, Rpg::FieldIo<3>, Cuda::RField<3> >;
} 
#endif

} // namespace Pscf
#endif
