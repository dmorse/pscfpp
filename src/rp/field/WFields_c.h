#ifndef RPC_W_FIELDS_H
#define RPC_W_FIELDS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2026, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/WFieldsBase.h>       // base class template
#include <pscf/backend/CPT.h>          // template argument
#include <pscf/backend/TmplDeclare.h>  // template declaration macros

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /**
   * A container of fields stored in both basis and r-grid format.
   *
   * Specializations of this template with D =1, 2, and 3 are derived
   * from specializations of the base class template Rp::WFieldsBase, 
   * and inherit their public interface and all of their source code
   * from this base class. This partial specialization is designed 
   * for use with a serial CPU backend (T=CPT).
   *
   * \ingroup Rp_Field_Module
   */
   template <int D>
   class WFields<D,CPT> 
    : public Rp::WFieldsBase<D,CPT>
   {};

   // Explicit instantiation declarations
   PSCF_TMPL_DECLARE_CPP(WFields);

} // namespace Rp
} // namespace Pscf
#endif
