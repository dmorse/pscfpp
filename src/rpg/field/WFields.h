#ifndef RPG_W_FIELDS_H
#define RPG_W_FIELDS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/WFields.h>    // base class template
#include <prdc/cuda/RField.h>    // base class member

namespace Pscf {
namespace Rpg {

   // Forward declaration
   template <int D> class FieldIo;

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cuda;

   /**
   * A container of fields stored in both basis and r-grid format.
   *
   * Almost all of the implementation of this class is defined by the base
   * class template Prdc::WFieldsTmpl . See documentation of that base
   * class template for documentation of most member functions.
   *
   * \ingroup Rpg_Field_Module
   */
   template <int D>
   class WFields
    : public Rp::WFields<D, RField<D>, FieldIo<D> >
   {

   public:

      /**
      * Set new w fields, in unfolded real-space (r-grid) format.
      *
      * The array fields is an unfolded array that contains fields for
      * all monomer types, with the field for monomer 0 first, etc.
      *
      * \param fields  unfolded array of new w (chemical potential) fields
      */
      void setRGrid(DeviceArray<cudaReal>& fields);

      // Alias for base class.
      using Base = Rp::WFields<D, RField<D>, FieldIo<D> >;

      // Declaration to prevent base class method from being hidden
      using Base::setRGrid;

   };

} // namespace Rpg
} // namespace Pscf

// Explicit instantiation declarations
namespace Pscf {
   namespace Rp {
      using namespace Prdc::Cuda;
      extern template class WFields<1, RField<1>, Rpg::FieldIo<1> >;
      extern template class WFields<2, RField<2>, Rpg::FieldIo<2> >;
      extern template class WFields<3, RField<3>, Rpg::FieldIo<3> >;
   }
   namespace Rpg {
      extern template class WFields<1>;
      extern template class WFields<2>;
      extern template class WFields<3>;
   }
}
#endif
