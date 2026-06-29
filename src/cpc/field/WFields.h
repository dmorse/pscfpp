#ifndef CPC_W_FIELDS_H
#define CPC_W_FIELDS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <cp/field/WFields.h>          // base class template
#include <prdc/field/cpu/CField.h>     // base class template argument
#include <cpc/field/FieldIo.h>   // base class template argument

namespace Pscf {
namespace Cpc {

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * A container of w fields.
   *
   * The public interface of this class is identical to that of the base
   * class template Pscf::Cp::WFields. Please see documentation
   * of that base class for API documentation.
   *
   * \ingroup Cpc_Field_Module
   */
   template <int D>
   class WFields 
     : public Cp::WFields<D, CField<D>, FieldIo<D> >
   {
   public:

      /// Alias for base class.
      using Base = Cp::WFields< D, CField<D>, FieldIo<D> >;

      // Inherited public member functions
      using Base::setFieldIo;
      using Base::allocate;
      using Base::setFields;
      using Base::clear;
      using Base::fields;
      using Base::field;
      using Base::isAllocated;
      using Base::hasData;

   protected:

      using Base::meshDimensions;
      using Base::meshSize;
      using Base::nMonomer;
      using Base::fieldIo;

   private:

      /**
      * Assign one CField<D> to another: lhs = rhs.
      *
      * \param lhs  left-hand side of assignment
      * \param rhs  right-hand side of assignment
      */
      void assignField(CField<D>& lhs, CField<D> const & rhs) const 
      override;

   };

   // Explicit instantiation declarations
   extern template class WFields<1>;
   extern template class WFields<2>;
   extern template class WFields<3>;

} // namespace Cpc
} // namespace Pscf

// Explicit instantiation declarations for base class
namespace Pscf {
   namespace Cp {
      extern template 
      class WFields<1, Prdc::Cpu::CField<1>, Cpc::FieldIo<1> >;
      extern template 
      class WFields<2, Prdc::Cpu::CField<2>, Cpc::FieldIo<2> >;
      extern template 
      class WFields<3, Prdc::Cpu::CField<3>, Cpc::FieldIo<3> >;
   } 
} 
#endif
