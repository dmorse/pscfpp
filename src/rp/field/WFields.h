#ifndef RP_W_FIELDS_H
#define RP_W_FIELDS_H

/*
* PSCF - Polymer Self-Consistent Field
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <rp/field/WFieldsBase.h>    // base class template
#include <prdc/field/TmplDeclare.h>  // explicit declaration macros

namespace Pscf {
namespace Rp {

   using namespace Util;
   using namespace Prdc;

   /**
   * A container of fields stored in both basis and r-grid format.
   *
   * This primary template is used for the CPU specialization, for which
   * T = CppTp<D>. A partial specializations is defined for the CUDA 
   * backend.
   *
   * \ingroup Rp_Field_Module
   */
   template <int D, class T>
   class WFields : public Rp::WFieldsBase<D,T>
   {};

   // Declare explicit specializations for CPU backend
   PRDC_TMPL_DECLARE_CPP(WFields);

} // namespace Rp
} // namespace Pscf

#ifdef PSCF_CUDA
#include <rp/field/WFields_cu.h>
#endif

#endif
