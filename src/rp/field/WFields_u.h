#ifndef RPG_W_FIELDS_H
#define RPG_W_FIELDS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/WFieldsBase.h>  // base class template
#include <pscf/backends/CUT.h>      // base class template argument

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   // Declare primary template
   template <int D, class T> class WFields;

   /**
   * A container of fields stored in both basis and r-grid format.
   *
   * Specializations of this template with D=1, 2, and 3 are derived from
   * corresponding specializations of base class template Rp::Class, and
   * inherit their public interface and almost all of their source code
   * from this base class.
   *
   * \see Rp::WFields
   * \ingroup Rp_Field_Module
   */
   template <int D>
   class WFields<D,CUT> 
    : public WFieldsBase<D,CUT>
   {

   public:

      WFields() = default;
      virtual ~WFields() = default;

      /**
      * Set new w fields, in unfolded real-space (r-grid) format.
      *
      * The input array fields is an unfolded array that contains fields 
      * for all monomer types, with the field for monomer 0 first, etc.
      *
      * \param fields  unfolded array of new w fields (input)
      */
      void setRGrid(DeviceArray<cudaReal>& fields);

      /// Alias for base class.
      using RpWFields = Rp::WFieldsBase<D,CUT>;

      // Declaration to avoid hiding overloaded base class method 
      using RpWFields::setRGrid;

   };

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE_CUDA(WFieldsBase);

} // namespace Rp
} // namespace Pscf
#endif
